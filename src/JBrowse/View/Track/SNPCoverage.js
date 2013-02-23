define( "JBrowse/View/Track/SNPCoverage", ['dojo/_base/declare',
         'dojo/_base/array',
         'JBrowse/View/Track/Wiggle',
         'JBrowse/Util',
         'JBrowse/Model/NestedFrequencyTable',
         'JBrowse/View/Track/MismatchesMixin'
        ],
        function( declare, array, Wiggle, Util, NestedFrequencyTable, MismatchesMixin ) {

var dojof = Util.dojof;

// feature class for the features we make for the calculated coverage
// values
var CoverageFeature = Util.fastDeclare(
    {
        get: function(f) { return this[f]; },
        tags: function() { return [ 'start', 'end', 'score' ]; },
        score: 0,
        constructor: function( args ) {
            this.start = args.start;
            this.end = args.end;
            this.score = args.score;
        }
    });


return declare( [Wiggle, MismatchesMixin],
{
    constructor: function() {
        // force conf variables that are meaningless for this kind of track, and maybe harmful
        delete this.config.bicolor_pivot;
        delete this.config.scale;
        delete this.config.align;
    },

    _defaultConfig: function() {
        return Util.deepUpdate(
            dojo.clone( this.inherited(arguments) ),
            {
                min_score: 0,
                max_score: 100
            }
        );
    },

    getGlobalStats: function() {
        return {};
    },

    getFeatures: function( query, featureCallback, finishCallback, errorCallback ) {
        var thisB = this;
        var leftBase  = query.start;
        var rightBase = query.end;
        var scale = query.scale; // px/bp
        var widthBp = rightBase-leftBase;
        var widthPx = widthBp * ( query.scale || 1/query.basesPerSpan);

        var binWidth = Math.ceil( query.basesPerSpan ); // in bp

        var binNumber = function( bp ) {
            return Math.floor( (bp-leftBase) / binWidth );
        };

        // init coverage bins
        var maxBin = binNumber( rightBase );
        var coverageBins = new Array( maxBin+1 );
        for( var i = 0; i <= maxBin; i++ ) {
            coverageBins[i] = new NestedFrequencyTable();
        }
        var binOverlap = function( bp, isRightEnd ) {
            var binCoord  = (bp-leftBase) / binWidth;
            var binNumber = Math.floor( binCoord );

            // only calculate the overlap if this lies in this block
            if( binNumber >= 0 && binNumber <= maxBin ) {
                var overlap =
                    isRightEnd ? 1 - ( binCoord - binNumber )
                               : binCoord - binNumber;
                return {
                    bin: binNumber,
                    overlap: overlap // between 0 and 1: proportion of this bin that the feature overlaps
                };
            }
            // otherwise null, this feature goes outside the block
            else {
                return isRightEnd ? { bin: maxBin, overlap: 1 }
                                  : { bin: 0,      overlap: 1 };
            }
        };


        thisB.store.getFeatures(
            query,
            function( feature ) {

                // calculate total coverage
                var startBO = binOverlap( feature.get('start'), false );
                var endBO   = binOverlap( feature.get('end')-1  , true  );

                // increment start and end partial-overlap bins by proportion of overlap
                if( startBO.bin == endBO.bin ) {
                    coverageBins[startBO.bin].increment( 'reference', endBO.overlap + startBO.overlap - 1 );
                }
                else {
                    coverageBins[startBO.bin].increment( 'reference', startBO.overlap );
                    coverageBins[endBO.bin].increment(   'reference', endBO.overlap   );
                }

                // increment completely overlapped interior bins by 1
                for( var i = startBO.bin+1; i <= endBO.bin-1; i++ ) {
                    coverageBins[i].increment( 'reference', 1 );
                }


                // Calculate SNP coverage
                if( binWidth == 1 ) {

                    // mark each bin as having its snps counted
                    for( var i = startBO.bin; i <= endBO.bin; i++ ) {
                        coverageBins[i].snpsCounted = 1;
                    }

                    // parse the MD
                    var mdTag = feature.get('MD');
                    if( mdTag ) {
                        var SNPs = thisB._mdToMismatches(feature, mdTag);
                        // loops through mismatches and updates coverage variables accordingly.
                        for (var i = 0; i<SNPs.length; i++) {
                            var pos = binNumber( feature.get('start') + SNPs[i].start );
                            var bin = coverageBins[pos];
                            if( bin ) {
                                var strand = { '-1': '-', '1': '+' }[ ''+feature.get('strand') ] || 'unstranded';
                                // Note: we decrement 'reference' so that total of the score is the total coverage
                                bin.decrement('reference', 1/binWidth );
                                var base = SNPs[i].bases;
                                bin.getNested(base).increment(strand, 1/binWidth);
                            }
                        }
                    }

                }
            },
            function () {
                var makeFeatures = function() {
                    // make fake features from the coverage
                    for( var i = 0; i <= maxBin; i++ ) {
                        var bpOffset = leftBase+binWidth*i;
                        featureCallback( new CoverageFeature({
                            start: bpOffset,
                            end:   bpOffset+binWidth,
                            score: coverageBins[i]
                         }));
                    }
                    finishCallback();
                };

                // if we are zoomed to base level, try to fetch the
                // reference sequence for this region and record each
                // of the bases in the coverage bins
                if( binWidth == 1 ) {
                    var sequence;
                    thisB.browser.getStore( 'refseqs', function( refSeqStore ) {
                        if( refSeqStore ) {
                            refSeqStore.getFeatures( query,
                                                     function(f) {
                                                         sequence = f.get('seq');
                                                     },
                                                     function() {
                                                         if( sequence ) {
                                                             for( var base = leftBase; base <= rightBase; base++ ) {
                                                                 var bin = binNumber( base );
                                                                 coverageBins[bin].refBase = sequence[bin];
                                                             }
                                                         }
                                                         makeFeatures();
                                                     },
                                                     makeFeatures
                                                   );
                        } else {
                            makeFeatures();
                        }
                    });
                } else {
                    makeFeatures();
                }
            }
        );
    },

    /*
     * Draw a set of features on the canvas.
     * @private
     */
    _drawFeatures: function( scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale ) {
        var thisB = this;
        var context = canvas.getContext('2d');
        var canvasHeight = canvas.height;
        var toY = dojo.hitch( this, function( val ) {
           return canvasHeight * ( 1-dataScale.normalize.call(this, val) );
        });
        var originY = toY( dataScale.origin );

        // a canvas element below the histogram that will contain indicators of likely SNPs
        var snpCanvasHeight = parseInt( block.parentNode.style.height ) - canvas.height;
        var snpCanvas = dojo.create('canvas',
                                    {height: snpCanvasHeight,
                                     width: canvas.width,
                                     style: {
                                         cursor: 'default',
                                         width: "100%",
                                         height: snpCanvasHeight + "px"
                                     },
                                     innerHTML: 'Your web browser cannot display this type of track.',
                                     className: 'SNP-indicator-track'
                                 }, block);
        var snpContext = snpCanvas.getContext('2d');

        var negColor  = this.config.style.neg_color;
        var clipColor = this.config.style.clip_marker_color;
        var bgColor   = this.config.style.bg_color;
        var disableClipMarkers = this.config.disable_clip_markers;

        var drawRectangle = function(ID, yPos, height, fRect) {
            // draw the background color if we are configured to do so
            if( bgColor && yPos >= 0 ) {
                context.fillStyle = bgColor;
                context.fillRect( fRect.l, 0, fRect.w, canvasHeight );
            }

            if( yPos <= canvasHeight ) { // if the rectangle is visible at all
                context.fillStyle = thisB.colorForBase(ID);
                if( yPos <= originY ) {
                    // bar goes upward
                    context.fillRect( fRect.l, yPos, fRect.w, height);
                    if( !disableClipMarkers && yPos < 0 ) { // draw clip marker if necessary
                        context.fillStyle = clipColor || negColor;
                        context.fillRect( fRect.l, 0, fRect.w, 2 );
                    }
                }
                else {
                    // bar goes downward
                    context.fillRect( fRect.l, originY, fRect.w, height );
                    if( !disableClipMarkers && yPos >= canvasHeight ) { // draw clip marker if necessary
                        context.fillStyle = clipColor || thisB.colorForBase(ID);
                        context.fillRect( fRect.l, canvasHeight-3, fRect.w, 2 );
                    }
                }
            }
        };

        dojo.forEach( features, function(f,i) {
            var fRect = featureRects[i];
            var score = f.get('score');
            var totalHeight = score.total();

            // draw indicators of SNPs if base coverage is greater than 50% of total coverage
            score.forEach( function( count, category ) {
                if ( category != 'reference' && count > 0.5*totalHeight ) {
                    snpContext.beginPath();
                    snpContext.arc( fRect.l + 0.5*fRect.w,
                                    0.40*snpCanvas.height,
                                    0.20*snpCanvas.height,
                                    1.75 * Math.PI,
                                    1.25 * Math.PI,
                                    false);
                    snpContext.lineTo(fRect.l + 0.5*fRect.w, 0);
                    snpContext.closePath();
                    snpContext.fillStyle = thisB.colorForBase(category);
                    snpContext.fill();
                    snpContext.lineWidth = 1;
                    snpContext.strokeStyle = 'black';
                    snpContext.stroke();
                }
            });

            // Note: 'reference' is done first to ensure the grey part of the graph is on top
            drawRectangle( 'reference', toY(totalHeight), originY-toY( score.get('reference'))+1, fRect);
            totalHeight -= score.get('reference');

            score.forEach( function( count, category ) {
                if ( category != 'reference' ) {
                    drawRectangle( category, toY(totalHeight), originY-toY( count )+1, fRect);
                    totalHeight -= count;
                }
            });
        }, this );
    },

    /**
     * parse a SAM MD tag to find mismatching bases of the template versus the reference
     * @returns {Array[Object]} array of mismatches and their positions
     */
    _mdToMismatches: function( feature, mdstring ) {
        var mismatchRecords = [];
        var curr = { start: 0, bases: '' };
        var seq = feature.get('seq');
        var nextRecord = function() {
              mismatchRecords.push( curr );
              curr = { start: curr.start + curr.bases.length, bases: ''};
        };
        array.forEach( mdstring.match(/(\d+|\^[a-z]+|[a-z])/ig), function( token ) {
          if( token.match(/^\d/) ) { // matching bases
              curr.start += parseInt( token );
          }
          else if( token.match(/^\^/) ) { // insertion in the template
              var i = token.length-1;
              while( i-- ) {
                  curr.bases = '*';
                  nextRecord();
              }
          }
          else if( token.match(/^[a-z]/i) ) { // mismatch
              curr.bases = seq.substr( curr.start, token.length );
              nextRecord();
          }
        });
        return mismatchRecords;
    },


    /*
     * The following method is required to override the equivalent method in "WiggleBase.js"
     * It displays more complete data.
     */
    _showPixelValue: function( scoreDisplay, score ) {
        if( ! score )
            return false;

        function fmtNum( num ) {
            return parseFloat( num ).toPrecision(6).replace(/0+$/,'').replace(/\.$/,'');
        }

        if( score.snpsCounted ) {
            var total = score.total();
            var scoreSummary = '<table>';
            function pctString( count ) {
                return Math.round(count/total*100)+'%';
            }
            scoreSummary +=
                  '<tr class="ref"><td>'
                + (score.refBase ? score.refBase+'*' : 'Ref')
                + '</td><td class="count">'
                + fmtNum( score.get('reference') )
                + '</td><td class="pct">'
                + pctString( score.get('reference') )
                + '</td></tr>';

            score.forEach( function( count, category ) {
                if( category == 'reference' ) return;

                // if this count has more nested categories, do counts of those
                var subdistribution = '';
                if( count.forEach ) {
                    subdistribution = [];
                    count.forEach( function( count, category ) {
                        subdistribution.push( fmtNum(count) + ' '+category );
                    });
                    subdistribution = subdistribution.join(', ');
                    if( subdistribution )
                        subdistribution = '('+subdistribution+')';
                }

                category = { '*': 'del' }[category] || category;
                scoreSummary += '<tr><td>'+category + '</td><td class="count">' + fmtNum(count) + '</td><td class="pct">'
                                   +pctString(count)+'</td><td class="subdist">'+subdistribution + '</td></tr>';
            });
            scoreSummary += '<tr class="total"><td>Total</td><td class="count">'+fmtNum(total)+'</td><td class="pct">&nbsp;</td><td class="subdist">&nbsp;</td></tr>';
            scoreDisplay.innerHTML = scoreSummary+'</table>';
            return true;
        } else {
            scoreDisplay.innerHTML = '<table><tr><td>Total</td><td class="count">'+fmtNum(score)+'</td></tr></table>';
            return true;
        }
    }
});
});