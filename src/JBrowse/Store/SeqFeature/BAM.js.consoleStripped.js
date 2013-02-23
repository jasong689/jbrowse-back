require({cache:{
'JBrowse/Store/SeqFeature/BAM/Util':function(){
define( [ 'jszlib/inflate',
          'jszlib/arrayCopy',
          'JBrowse/Util'
        ],
        function( inflate, arrayCopy, Util ) {

var VirtualOffset = Util.fastDeclare({
    constructor: function(b, o) {
        this.block = b;
        this.offset = o;
    },
    toString: function() {
        return '' + this.block + ':' + this.offset;
    },
    cmp: function(b) {
        var a = this;
        return b.block - a.block || b.offset - a.offset;
    }
});

/**
 * @lends JBrowse.Store.SeqFeature.BAM.Util
 * Package of utility functions used in various places in the BAM code.
 */
var Utils = {

    readInt: function(ba, offset) {
        return (ba[offset + 3] << 24) | (ba[offset + 2] << 16) | (ba[offset + 1] << 8) | (ba[offset]);
    },

    readShort: function(ba, offset) {
        return (ba[offset + 1] << 8) | (ba[offset]);
    },

    readFloat: function(ba, offset) {
        var temp = new Uint8Array( 4 );
        for( var i = 0; i<4; i++ ) {
            temp[i] = ba[offset+i];
        }
        var fa = new Float32Array( temp.buffer );
        return fa[0];
    },

    readVirtualOffset: function(ba, offset) {
        //0 && console.log( 'readVob', offset );
        var block = (ba[offset+6] & 0xff) * 0x100000000
            + (ba[offset+5] & 0xff) * 0x1000000
            + (ba[offset+4] & 0xff) * 0x10000
            + (ba[offset+3] & 0xff) * 0x100
            + (ba[offset+2] & 0xff);
        var bint = (ba[offset+1] << 8) | ba[offset];
        if (block == 0 && bint == 0) {
            return null;  // Should only happen in the linear index?
        } else {
            return new VirtualOffset(block, bint);
        }
    },

    unbgzf: function(data, lim) {
        lim = Math.min( lim || Infinity, data.byteLength - 27);
        var oBlockList = [];
        var totalSize = 0;

        for( var ptr = [0]; ptr[0] < lim; ptr[0] += 8) {

            var ba = new Uint8Array( data, ptr[0], 18 );

            // check the bgzf block magic
            if( !( ba[0] == 31 && ba[1] == 139 ) ) {
                console.error( 'invalid BGZF block header, skipping', ba );
                break;
            }

            var xlen = Utils.readShort( ba, 10 );
            var compressedDataOffset = ptr[0] + 12 + xlen;

            // var inPtr = ptr[0];
            // var bSize = Utils.readShort( ba, 16 );
            // var logLength = Math.min(data.byteLength-ptr[0], 40);
            // 0 && console.log( xlen, bSize, bSize - xlen - 19, new Uint8Array( data, ptr[0], logLength ), logLength );

            var unc;
            try {
                unc = inflate(
                    data,
                    compressedDataOffset,
                    data.byteLength - compressedDataOffset,
                    ptr
                );
            } catch( inflateError ) {
                // if we have a buffer error and we have already
                // inflated some data, there is probably just an
                // incomplete BGZF block at the end of the data, so
                // ignore it and stop inflating
                if( /^Z_BUF_ERROR/.test(inflateError.statusString) && oBlockList.length ) {
                    break;
                }
                // otherwise it's some other kind of real error
                else {
                    throw inflateError;
                }
            }
            if( unc.byteLength ) {
                totalSize += unc.byteLength;
                oBlockList.push( unc );
            }
            // else {
            //     console.error( 'BGZF decompression failed for block ', compressedDataOffset, data.byteLength-compressedDataOffset, [inPtr] );
            // }
        }

        if (oBlockList.length == 1) {
            return oBlockList[0];
        } else {
            var out = new Uint8Array(totalSize);
            var cursor = 0;
            for (var i = 0; i < oBlockList.length; ++i) {
                var b = new Uint8Array(oBlockList[i]);
                arrayCopy(b, 0, out, cursor, b.length);
                cursor += b.length;
            }
            return out.buffer;
        }
    }
};

return Utils;

});
},
'JBrowse/Store/SeqFeature/BAM/File':function(){
define( [
            'dojo/_base/declare',
            'dojo/_base/array',
            'JBrowse/Util',
            'JBrowse/Store/LRUCache',
            './Util',
            './LazyFeature'
        ],
        function( declare, array, Util, LRUCache, BAMUtil, BAMFeature ) {

var BAM_MAGIC = 21840194;
var BAI_MAGIC = 21578050;

var dlog = function(){ console.error.apply(console, arguments); };

//
// Binning (transliterated from SAM1.3 spec)
//

/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
function reg2bin(beg, end)
{
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
    return 0;
}

/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
var MAX_BIN = (((1<<18)-1)/7);
function reg2bins(beg, end)
{
    var k, list = [];
    --end;
    list.push(0);
    for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) list.push(k);
    for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) list.push(k);
    for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) list.push(k);
    for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) list.push(k);
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list.push(k);
    return list;
}

var Chunk = Util.fastDeclare({
    constructor: function(minv,maxv,bin) {
        this.minv = minv;
        this.maxv = maxv;
        this.bin = bin;
    },
    toString: function() {
        return this.minv+'..'+this.maxv+' (bin '+this.bin+')';
    }
});

var readInt   = BAMUtil.readInt;
var readVirtualOffset = BAMUtil.readVirtualOffset;

var BamFile = declare( null,


/**
 * @lends JBrowse.Store.SeqFeature.BAM.File
 */
{

    /**
     * Low-level BAM file reading code.
     *
     * Adapted by Robert Buels from bam.js in the Dalliance Genome
     * Explorer which is copyright Thomas Down 2006-2010
     * @constructs
     */
    constructor: function( args ) {
        this.store = args.store;
        this.data  = args.data;
        this.bai   = args.bai;
    },

    init: function( args ) {
        var bam = this;
        var successCallback = args.success || function() {};
        var failCallback = args.failure || function() {};

        this._readBAI( dojo.hitch( this, function() {
            this._readBAMheader( function() {
                successCallback();
            }, failCallback );
        }), failCallback );
    },

    _readBAI: function( successCallback, failCallback ) {
        // Do we really need to fetch the whole thing? :-(
        this.bai.fetch( dojo.hitch( this, function(header) {
            if (!header) {
                dlog("No data read from BAM index (BAI) file");
                failCallback("No data read from BAM index (BAI) file");
                return;
            }

            if( ! Uint8Array ) {
                dlog('Browser does not support typed arrays');
                failCallback('Browser does not support typed arrays');
                return;
            }

            var uncba = new Uint8Array(header);
            if( readInt(uncba, 0) != BAI_MAGIC) {
                dlog('Not a BAI file');
                failCallback('Not a BAI file');
                return;
            }

            var nref = readInt(uncba, 4);

            this.indices = [];

            var p = 8;
            for (var ref = 0; ref < nref; ++ref) {
                var blockStart = p;
                var nbin = readInt(uncba, p); p += 4;
                for (var b = 0; b < nbin; ++b) {
                    var bin = readInt(uncba, p);
                    var nchnk = readInt(uncba, p+4);
                    p += 8;
                    for( var chunkNum = 0; chunkNum < nchnk; chunkNum++ ) {
                        var vo = readVirtualOffset( uncba, p );
                        this._findMinAlignment( vo );
                        p += 16;
                    }
                }
                var nintv = readInt(uncba, p); p += 4;
                // as we're going through the linear index, figure out
                // the smallest virtual offset in the indexes, which
                // tells us where the BAM header ends
                this._findMinAlignment( nintv ? readVirtualOffset(uncba,p) : null );

                p += nintv * 8;
                if( nbin > 0 || nintv > 0 ) {
                    this.indices[ref] = new Uint8Array(header, blockStart, p - blockStart);
                }
            }

            this.empty = ! this.indices.length;

            successCallback( this.indices, this.minAlignmentVO );
        }), failCallback );
    },

    _findMinAlignment: function( candidate ) {
        if( candidate && ( ! this.minAlignmentVO || this.minAlignmentVO.cmp( candidate ) < 0 ) )
            this.minAlignmentVO = candidate;
    },

    _readBAMheader: function( successCallback, failCallback ) {
        var thisB = this;
        // We have the virtual offset of the first alignment
        // in the file.  Cannot completely determine how
        // much of the first part of the file to fetch to get just
        // up to that, since the file is compressed.  Thus, fetch
        // up to the start of the BGZF block that the first
        // alignment is in, plus 64KB, which should get us that whole
        // BGZF block, assuming BGZF blocks are no bigger than 64KB.
        thisB.data.read(
            0,
            thisB.minAlignmentVO ? thisB.minAlignmentVO.block + 65535 : null,
            function(r) {
                var unc = BAMUtil.unbgzf(r);
                var uncba = new Uint8Array(unc);

                if( readInt(uncba, 0) != BAM_MAGIC) {
                    dlog('Not a BAM file');
                    failCallback( 'Not a BAM file' );
                    return;
                }

                var headLen = readInt(uncba, 4);
                // var header = '';
                // for (var i = 0; i < headLen; ++i) {
                //     header += String.fromCharCode(uncba[i + 8]);
                // }


                // have to do another request, because sometimes
                // minAlignment VO is just flat wrong.
                // if headLen is not too big, this will just be in the
                // RemoteBinaryFile cache
                thisB.data.read( 0, headLen+8+65536,
                                 function(r) {
                    var unc = BAMUtil.unbgzf(r);
                    var uncba = new Uint8Array(unc);

                    var nRef = readInt(uncba, headLen + 8);
                    var p = headLen + 12;

                    thisB.chrToIndex = {};
                    thisB.indexToChr = [];
                    for (var i = 0; i < nRef; ++i) {
                        var lName = readInt(uncba, p);
                        var name = '';
                        for (var j = 0; j < lName-1; ++j) {
                            name += String.fromCharCode(uncba[p + 4 + j]);
                        }
                        var lRef = readInt(uncba, p + lName + 4);
                        // dlog(name + ': ' + lRef);
                        thisB.chrToIndex[name] = i;
                        if (name.indexOf('chr') == 0) {
                            thisB.chrToIndex[name.substring(3)] = i;
                        } else {
                            thisB.chrToIndex['chr' + name] = i;
                        }

                        thisB.indexToChr.push({ name: name, length: lRef });

                        p = p + 8 + lName;
                    }

                    successCallback();
            }, failCallback );
        }, failCallback );
    },

    blocksForRange: function(refId, min, max) {
        var index = this.indices[refId];
        if (!index) {
            return [];
        }

        var intBinsL = reg2bins(min, max);
        var intBins = [];
        for (var i = 0; i < intBinsL.length; ++i) {
            intBins[intBinsL[i]] = true;
        }
        var leafChunks = [], otherChunks = [];

        var nbin = readInt(index, 0);
        var p = 4;
        for (var b = 0; b < nbin; ++b) {
            var bin = readInt(index, p);
            var nchnk = readInt(index, p+4);
    //        dlog('bin=' + bin + '; nchnk=' + nchnk);
            p += 8;
            if (intBins[bin]) {
                for (var c = 0; c < nchnk; ++c) {
                    var cs = readVirtualOffset(index, p);
                    var ce = readVirtualOffset(index, p + 8);
                    (bin < 4681 ? otherChunks : leafChunks).push(new Chunk(cs, ce, bin));
                    p += 16;
                }
            } else {
                p +=  (nchnk * 16);
            }
        }
    //    dlog('leafChunks = ' + miniJSONify(leafChunks));
    //    dlog('otherChunks = ' + miniJSONify(otherChunks));

        var nintv = readInt(index, p);
        var lowest = null;
        var minLin = Math.min(min>>14, nintv - 1), maxLin = Math.min(max>>14, nintv - 1);
        for (var i = minLin; i <= maxLin; ++i) {
            var lb =  readVirtualOffset(index, p + 4 + (i * 8));
            if (!lb) {
                continue;
            }
            if (!lowest || lb.block < lowest.block || lb.offset < lowest.offset) {
                lowest = lb;
            }
        }
        // dlog('Lowest LB = ' + lowest);

        var prunedOtherChunks = [];
        if (lowest != null) {
            for (var i = 0; i < otherChunks.length; ++i) {
                var chnk = otherChunks[i];
                if (chnk.maxv.block >= lowest.block && chnk.maxv.offset >= lowest.offset) {
                    prunedOtherChunks.push(chnk);
                }
            }
        }
        // dlog('prunedOtherChunks = ' + miniJSONify(prunedOtherChunks));
        otherChunks = prunedOtherChunks;

        var intChunks = [];
        for (var i = 0; i < otherChunks.length; ++i) {
            intChunks.push(otherChunks[i]);
        }
        for (var i = 0; i < leafChunks.length; ++i) {
            intChunks.push(leafChunks[i]);
        }

        intChunks.sort(function(c0, c1) {
            var dif = c0.minv.block - c1.minv.block;
            if (dif != 0) {
                return dif;
            } else {
                return c0.minv.offset - c1.minv.offset;
            }
        });
        var mergedChunks = [];
        if (intChunks.length > 0) {
            var cur = intChunks[0];
            for (var i = 1; i < intChunks.length; ++i) {
                var nc = intChunks[i];
                if (nc.minv.block == cur.maxv.block /* && nc.minv.offset == cur.maxv.offset */) { // no point splitting mid-block
                    cur = new Chunk(cur.minv, nc.maxv, 'linear');
                } else {
                    mergedChunks.push(cur);
                    cur = nc;
                }
            }
            mergedChunks.push(cur);
        }
    //    dlog('mergedChunks = ' + miniJSONify(mergedChunks));

        return mergedChunks;
    },

    fetch: function(chr, min, max, callback) {

        var chrId = this.chrToIndex && this.chrToIndex[chr];
        var chunks;
        if( !( chrId >= 0 ) ) {
            chunks = [];
        } else {
            chunks = this.blocksForRange(chrId, min, max);
            if (!chunks) {
                callback(null, 'Error in index fetch');
            }
        }

        // toString function is used by the cache for making cache keys
        chunks.toString = function() {
            return this.join(', ');
        };

        try {
            this._fetchChunkFeatures( chunks, function( features, error ) {
                if( error ) {
                    callback( null, error );
                } else {
                    features = array.filter( features, function( feature ) {
                        return ( !( feature.get('end') < min || feature.get('start') > max )
                                 && ( chrId === undefined || feature._refID == chrId ) );
                    });
                    callback( features );
                }
            });
        } catch( e ) {
            callback( null, e );
        }
    },

    _fetchChunkFeatures: function( chunks, callback ) {
        var thisB = this;

        if( ! chunks.length ) {
            callback([]);
            return;
        }

        var features = [];
        var chunksProcessed = 0;

        var cache = this.featureCache = this.featureCache || new LRUCache({
            name: 'bamFeatureCache',
            fillCallback: dojo.hitch( this, '_readChunk' ),
            sizeFunction: function( features ) {
                return features.length;
            },
            maxSize: 100000 // cache up to 100,000 BAM features
        });

        var error;
        array.forEach( chunks, function( c ) {
            cache.get( c, function( f, e ) {
                error = error || e;
                features.push.apply( features, f );
                if( ++chunksProcessed == chunks.length )
                    callback( features, error );
            });
        });

    },

    _readChunk: function( chunk, callback ) {
        var thisB = this;
        var features = [];
        var fetchMin = chunk.minv.block;
        var fetchMax = chunk.maxv.block + (1<<16); // *sigh*

        thisB.data.read(fetchMin, fetchMax - fetchMin + 1, function(r) {
            try {
                var data = BAMUtil.unbgzf(r, chunk.maxv.block - chunk.minv.block + 1);
                thisB.readBamFeatures( new Uint8Array(data), chunk.minv.offset, features, callback );
            } catch( e ) {
                callback( null, e );
            }
        });
    },

    readBamFeatures: function(ba, blockStart, sink, callback ) {
        var that = this;
        var featureCount = 0;

        var maxFeaturesWithoutYielding = 300;

        while ( true ) {
            if( blockStart >= ba.length ) {
                // if we're done, call the callback and return
                callback( sink );
                return;
            }
            else if( featureCount <= maxFeaturesWithoutYielding ) {
                // if we've read no more than 200 features this cycle, read another one
                var blockSize = readInt(ba, blockStart);
                var blockEnd = blockStart + 4 + blockSize - 1;

                // only try to read the feature if we have all the bytes for it
                if( blockEnd < ba.length ) {
                    var feature = new BAMFeature({
                        store: this.store,
                        file: this,
                        bytes: { byteArray: ba, start: blockStart, end: blockEnd }
                     });
                    sink.push(feature);
                    featureCount++;
                }

                blockStart = blockEnd+1;
            }
            else {
                // if we're not done but we've read a good chunk of
                // features, put the rest of our work into a timeout to continue
                // later, avoiding blocking any UI stuff that's going on
                window.setTimeout( function() {
                    that.readBamFeatures( ba, blockStart, sink, callback );
                }, 1);
                return;
            }
        }
    }
});

return BamFile;

});

},
'JBrowse/Store/SeqFeature/BAM/LazyFeature':function(){
define( ['dojo/_base/array',
         'JBrowse/Util',
         './Util',
         'JBrowse/Model/SimpleFeature'
        ],
        function( array, Util, BAMUtil, SimpleFeature ) {

var SEQRET_DECODER = ['=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'];
var CIGAR_DECODER  = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'];

var readInt   = BAMUtil.readInt;
var readShort = BAMUtil.readShort;
var readFloat = BAMUtil.readFloat;

var Feature = Util.fastDeclare(
{
    constructor: function( args ) {
        this.store = args.store;
        this.file  = args.file;
        this.data  = {
            type: 'match',
            source: args.store.source
        };
        this.bytes = {
            start: args.bytes.start,
            end: args.bytes.end,
            byteArray: args.bytes.byteArray
        };

        this._coreParse();
    },

    get: function( field) {
        return this._get( field.toLowerCase() );
    },

    // same as get(), except requires lower-case arguments.  used
    // internally to save lots of calls to field.toLowerCase()
    _get: function( field ) {
        return field in this.data ? this.data[field] : // have we already parsed it out?
            function() {
                var v = this.data[field] =
                    this[field]            ? this[field]()            : // maybe we have a special parser for it
                    this._flagMasks[field] ? this._parseFlag( field ) : // or is it a flag?
                                             this._parseTag( field );   // otherwise, look for it in the tags
                return v;
            }.call(this);
    },

    tags: function() {
        return this._get('_tags');
    },

    _tags: function() {
        this._parseAllTags();

        var tags = [ 'seq', 'seq_reverse_complemented', 'unmapped' ];
        if( ! this._get('unmapped') )
            tags.push( 'start', 'end', 'strand', 'score', 'qual', 'MQ', 'CIGAR', 'length_on_ref' );
        if( this._get('multi_segment_template') ) {
            tags.push( 'multi_segment_all_aligned',
                       'multi_segment_next_segment_unmapped',
                       'multi_segment_next_segment_reversed',
                       'multi_segment_first',
                       'multi_segment_last',
                       'secondary_alignment',
                       'qc_failed',
                       'duplicate',
                       'next_segment_position'
                     );
        }
        tags = tags.concat( this._tagList || [] );

        var d = this.data;
        for( var k in d ) {
            if( d.hasOwnProperty( k ) && k[0] != '_' )
                tags.push( k );
        }

        var seen = {};
        tags = array.filter( tags, function(t) {
            if( t in this.data && this.data[t] === undefined )
                return false;

            var lt = t.toLowerCase();
            var s = seen[lt];
            seen[lt] = true;
            return ! s;
        },this);

        return tags;
    },

    parent: function() {
        return undefined;
    },

    children: function() {
        return this._get('subfeatures');
    },

    id: function() {
        return this._get('name')+'/'+this._get('md')+'/'+this._get('cigar')+'/'+this._get('start');
    },

    // special parsers
    /**
     * Mapping quality score.
     */
    mq: function() {
        var mq = (this._get('_bin_mq_nl') & 0xff00) >> 8;
        return mq == 255 ? undefined : mq;
    },
    score: function() {
        return this._get('mq');
    },
    qual: function() {
        if( this._get('unmapped') )
            return undefined;

        var qseq = [];
        var byteArray = this.bytes.byteArray;
        var p = this.bytes.start + 36 + this._get('_l_read_name') + this._get('_n_cigar_op')*4 + this._get('_seq_bytes');
        var lseq = this._get('seq_length');
        for (var j = 0; j < lseq; ++j) {
            qseq.push( byteArray[p + j] );
        }
        return qseq.join(' ');
    },
    strand: function() {
        var xs = this._get('xs');
        return xs ? ( xs == '-' ? -1 : 1 ) :
               this._get('seq_reverse_complemented') ? -1 :  1;
    },
    /**
     * Length in characters of the read name.
     */
    _l_read_name: function() {
        return this._get('_bin_mq_nl') & 0xff;
    },
    /**
     * number of bytes in the sequence field
     */
    _seq_bytes: function() {
        return (this._get('seq_length') + 1) >> 1;
    },
    seq: function() {
        var seq = '';
        var byteArray = this.bytes.byteArray;
        var p = this.bytes.start + 36 + this._get('_l_read_name') + this._get('_n_cigar_op')*4;
        var seqBytes = this._get('_seq_bytes');
        for (var j = 0; j < seqBytes; ++j) {
            var sb = byteArray[p + j];
            seq += SEQRET_DECODER[(sb & 0xf0) >> 4];
            seq += SEQRET_DECODER[(sb & 0x0f)];
        }
        return seq;
    },
    name: function() {
        return this._get('_read_name');
    },
    _read_name: function() {
        var byteArray = this.bytes.byteArray;
        var readName = '';
        var nl = this._get('_l_read_name');
        var p = this.bytes.start + 36;
        for (var j = 0; j < nl-1; ++j) {
            readName += String.fromCharCode(byteArray[p+j]);
        }
        return readName;
    },
    _n_cigar_op: function() {
        return this._get('_flag_nc') & 0xffff;
    },
    cigar: function() {
        if( this._get('unmapped') )
            return undefined;

        var byteArray   = this.bytes.byteArray;
        var numCigarOps = this._get('_n_cigar_op');
        var p = this.bytes.start + 36 + this._get('_l_read_name');
        var cigar = '';
        var lref = 0;
        for (var c = 0; c < numCigarOps; ++c) {
            var cigop = readInt(byteArray, p);
            var lop = cigop >> 4;
            var op = CIGAR_DECODER[cigop & 0xf];
            cigar = cigar + lop + op;
            lref += lop;
            p += 4;
        }
        this.data.length_on_ref = lref;
        return cigar;
    },
    next_segment_position: function() {
        var nextRefID = this._get('_next_refid');
        var nextSegment = this.file.indexToChr[nextRefID];
        if( nextSegment )
            return nextSegment.name+':'+this._get('_next_pos');
        else
            return undefined;
    },
    subfeatures: function() {
        if( ! this.store.createSubfeatures )
            return undefined;

        var cigar = this._get('cigar');
        if( cigar )
            return this._cigarToSubfeats( cigar );

        return undefined;
    },
    length_on_ref: function() {
        var c = this._get('cigar'); // the length_on_ref is set as a
                                   // side effect of the CIGAR parsing
        return this.data.length_on_ref;
    },
    _flags: function() {
        return (this.data._flag_nc & 0xffff0000) >> 16;
    },
    end: function() {
        return this._get('start') + ( this._get('length_on_ref') || this._get('seq_length') || undefined );
    },

    seq_id: function() {
        if( this._get('unmapped') )
            return undefined;

        return ( this.file.indexToChr[ this._refID ] || {} ).name;
    },

    /**
     * parse the core data: start, end, strand, etc
     */
    _coreParse: function() {
        var byteArray = this.bytes.byteArray;
        var dataStart = this.bytes.start+4;

        var tempBytes = new Uint8Array( 32 );
        for( var i = 0; i<32; i++ ) {
            tempBytes[i] = byteArray[i+dataStart];
        }
        var ints = new Int32Array( tempBytes.buffer );

        var d = this.data;
        this._refID        = ints[0];
        d.start            = ints[1];
        d._bin_mq_nl       = ints[2];
        d._flag_nc         = ints[3];
        d.seq_length       = ints[4];
        d._next_refid      = ints[5];
        d._next_pos        = ints[6];
        d.template_length  = ints[7];
    },

    /**
     * Get the value of a tag, parsing the tags as far as necessary.
     * Only called if we have not already parsed that field.
     */
    _parseTag: function( tagName ) {
        // if all of the tags have been parsed and we're still being
        // called, we already know that we have no such tag, because
        // it would already have been cached.
        if( this._allTagsParsed )
            return undefined;

        this._tagList = this._tagList || [];
        var byteArray = this.bytes.byteArray;
        var p = this._tagOffset || this.bytes.start + 36 + this._get('_l_read_name') + this._get('_n_cigar_op')*4 + this._get('_seq_bytes') + this._get('seq_length');

        var blockEnd = this.bytes.end;
        while( p < blockEnd && lcTag != tagName ) {
            var tag      = String.fromCharCode( byteArray[p], byteArray[ p+1 ] );
            var lcTag    = tag.toLowerCase();
            var type = String.fromCharCode( byteArray[ p+2 ] );
            p += 3;

            var value;
            switch( type.toLowerCase() ) {
            case 'a':
                value = String.fromCharCode( byteArray[p] );
                p += 1;
                break;
            case 'i':
                value = readInt(byteArray, p );
                p += 4;
                break;
            case 'c':
                value = byteArray[p];
                p += 1;
                break;
            case 's':
                value = readShort(byteArray, p);
                p += 2;
                break;
            case 'f':
                value = readFloat( byteArray, p );
                p += 4;
                break;
            case 'z':
            case 'h':
                value = '';
                while( p <= blockEnd ) {
                    var cc = byteArray[p++];
                    if( cc == 0 ) {
                        break;
                    }
                    else {
                        value += String.fromCharCode(cc);
                    }
                }
                break;
            default:
                console.warn( "Unknown BAM tag type '"+type
                              +"', tags may be incomplete"
                            );
                value = undefined;
                p = blockEnd; // stop parsing tags
            }

            this._tagOffset = p;

            this._tagList.push( tag );
            if( lcTag == tagName )
                return value;
            else {
                this.data[ lcTag ] = value;
            }
        }
        this._allTagsParsed = true;
        return undefined;
    },
    _parseAllTags: function() {
        this._parseTag(); // calling _parseTag with no arg just parses
        // all the tags and returns the last one
    },

    _flagMasks: {
        multi_segment_template:              0x1,
        multi_segment_all_aligned:           0x2,
        unmapped:                            0x4,
        multi_segment_next_segment_unmapped: 0x8,
        seq_reverse_complemented:            0x10,
        multi_segment_next_segment_reversed: 0x20,
        multi_segment_first:                 0x40,
        multi_segment_last:                  0x80,
        secondary_alignment:                 0x100,
        qc_failed:                           0x200,
        duplicate:                           0x400
    },

    _parseFlag: function( flagName ) {
        return !!( this._get('_flags') & this._flagMasks[flagName] );
    },

    _parseCigar: function( cigar ) {
        return array.map( cigar.match(/\d+\D/g), function( op ) {
           return [ op.match(/\D/)[0].toUpperCase(), parseInt( op ) ];
        });
    },

    /**
     *  take a cigar string, and initial position, return an array of subfeatures
     */
    _cigarToSubfeats: function(cigar)    {
        var subfeats = [];
        var min = this._get('start');
        var max;
        var ops = this._parseCigar( cigar );
        for (var i = 0; i < ops.length; i++)  {
            var lop = ops[i][1];
            var op = ops[i][0];  // operation type
            // converting "=" to "E" to avoid possible problems later with non-alphanumeric type name
            if (op === "=")  { op = "E"; }

            switch (op) {
            case 'M':
            case 'D':
            case 'N':
            case 'E':
            case 'X':
                max = min + lop;
                break;
            case 'I':
                max = min;
                break;
            case 'P':  // not showing padding deletions (possibly change this later -- could treat same as 'I' ?? )
            case 'H':  // not showing hard clipping (since it's unaligned, and offset arg meant to be beginning of aligned part)
            case 'S':  // not showing soft clipping (since it's unaligned, and offset arg meant to be beginning of aligned part)
                break;
                // other possible cases
            }
            if( op !== 'N' ) {
                var subfeat = new SimpleFeature({
                    data: {
                    type: op,
                        start: min,
                        end: max,
                        strand: this._get('strand'),
                        cigar_op: lop+op
                    },
                    parent: this
                });
                subfeats.push(subfeat);
            }
            min = max;
        }
        return subfeats;
    }

});

return Feature;
});
},
'JBrowse/Model/SimpleFeature':function(){
define([
        'JBrowse/Util'
       ],
       function( Util ) {

var counter = 0;

return Util.fastDeclare({

    constructor: function( args ) {
        args = args || {};
        this.data = args.data || {};
        this._uniqueID = args.id || 'SimpleFeature_'+(counter++);
        this._parent = args.parent;
    },

    get: function(name) {
        return this.data[ name ];
    },

    set: function( name, val ) {
        this.data[ name ] = val;
    },

    tags: function() {
        var t = [];
        var d = this.data;
        for( var k in d ) {
            if( d.hasOwnProperty( k ) )
                t.push( k );
        }
        return t;
    },

    id: function( newid ) {
        if( newid )
            this._uniqueID = newid;
        return this._uniqueID;
    },

    parent: function() {
        return this._parent;
    },

    children: function() {
        return this.get('subfeatures');
    }

});

});
}}});
define( "JBrowse/Store/SeqFeature/BAM", [
            'dojo/_base/declare',
            'dojo/_base/array',
            'dojo/_base/Deferred',
            'dojo/_base/lang',
            'JBrowse/Util',
            'JBrowse/Store/SeqFeature',
            'JBrowse/Store/DeferredStatsMixin',
            'JBrowse/Store/DeferredFeaturesMixin',
            'JBrowse/Model/XHRBlob',
            './BAM/Util',
            './BAM/File'
        ],
        function( declare, array, Deferred, lang, Util, SeqFeatureStore, DeferredStatsMixin, DeferredFeaturesMixin, XHRBlob, BAMUtil, BAMFile ) {

var BAMStore = declare( [ SeqFeatureStore, DeferredStatsMixin, DeferredFeaturesMixin ],

/**
 * @lends JBrowse.Store.SeqFeature.BAM
 */
{
    /**
     * Data backend for reading feature data directly from a
     * web-accessible BAM file.
     *
     * @constructs
     */
    constructor: function( args ) {

        this.createSubfeatures = args.subfeatures;

        var bamBlob = args.bam || (function() {
                                       var url = Util.resolveUrl(
                                           args.baseUrl || '/',
                                           Util.fillTemplate( args.urlTemplate || 'data.bam',
                                           {'refseq': (this.refSeq||{}).name }
                                                            )
                                       );
                                       return new XHRBlob( url );
                                   }).call(this);
        var baiBlob = args.bai || (function() {
                                      var url = Util.resolveUrl(
                                          args.baseUrl || '/',
                                          Util.fillTemplate( args.baiUrlTemplate || args.urlTemplate+'.bai' || 'data.bam.bai',
                                                             {'refseq': (this.refSeq||{}).name }
                                                           )
                                      );
                                      return new XHRBlob( url );
                                  }).call(this);

        this.bam = new BAMFile({
                store: this,
                data: bamBlob,
                bai: baiBlob
        });

        this.source = ( bamBlob.url  ? bamBlob.url.match( /\/([^/\#\?]+)($|[\#\?])/ )[1] :
                        bamBlob.blob ? bamBlob.blob.name : undefined ) || undefined;

        this.bam.init({
            success: dojo.hitch( this, '_estimateGlobalStats',
                                 dojo.hitch( this, function(error) {
                                     if( error )
                                         this._failAllDeferred( error );
                                     else {
                                         this._deferred.stats.resolve({success:true});
                                         this._deferred.features.resolve({success:true});
                                     }

                                 })),
            failure: dojo.hitch( this, '_failAllDeferred' )
        });
    },


    /**
     * Fetch a region of the current reference sequence and use it to
     * estimate the feature density in the BAM file.
     * @private
     */
    _estimateGlobalStats: function( finishCallback ) {

        var statsFromInterval = function( refSeq, length, callback ) {
            var sampleCenter = refSeq.start*0.75 + refSeq.end*0.25;
            var start = Math.max( 0, Math.round( sampleCenter - length/2 ) );
            var end = Math.min( Math.round( sampleCenter + length/2 ), refSeq.end );
            this.bam.fetch( refSeq.name, start, end, dojo.hitch( this, function( features, error) {
                if ( error ) {
                    console.error( error );
                    callback.call( this, length,  null, error );
                }
                else if( features ) {
                    features = array.filter( features, function(f) { return f.get('start') >= start && f.get('end') <= end; } );
                    callback.call( this, length,
                                   {
                                       featureDensity: features.length / length,
                                       _statsSampleFeatures: features.length,
                                       _statsSampleInterval: { ref: refSeq.name, start: start, end: end, length: length }
                                   });
                }
            }));
        };

        var maybeRecordStats = function( interval, stats, error ) {
            if( error ) {
                finishCallback( error );
            } else {
                var refLen = this.refSeq.end - this.refSeq.start;
                 if( stats._statsSampleFeatures >= 300 || interval * 2 > refLen || error ) {
                     this.globalStats = stats;
                     0 && console.log( 'BAM statistics: '+this.source, stats );
                     finishCallback();
                 } else {
                     statsFromInterval.call( this, this.refSeq, interval * 2, maybeRecordStats );
                 }
            }
        };

        statsFromInterval.call( this, this.refSeq, 100, maybeRecordStats );
    },


    // called by getFeatures from the DeferredFeaturesMixin
    _getFeatures: function( query, featCallback, endCallback, errorCallback ) {
        var start = query.start;
        var end   = query.end;

        var maxFeaturesWithoutYielding = 300;
        this.bam.fetch( this.refSeq.name, start, end, function( features, error) {
                if ( error ) {
                    console.error( 'error fetching BAM data: ' + error );
                    if( errorCallback ) errorCallback( error );
                    return;
                }
                if( features ) {
                    var i = 0;
                    var readFeatures = function() {
                        for( ; i < features.length; i++ ) {
                            var feature = features[i];
                            // skip if this alignment is unmapped, or if it does not actually overlap this range
                            if (! (feature.get('unmapped') || feature.get('end') <= start || feature.get('start') >= end) )
                                try {
                                    featCallback( feature );
                                } catch(e) {
                                    if( errorCallback )
                                        errorCallback( e );
                                    else
                                        console.error( e, e.stack );
                                    return;
                                }

                            if( i && !( i % maxFeaturesWithoutYielding ) ) {
                                window.setTimeout( readFeatures, 1 );
                                i++;
                                break;
                            }
                        }
                        if( i >= features.length )
                            endCallback();
                    };

                    readFeatures();

                }
            });
    }

});

return BAMStore;
});