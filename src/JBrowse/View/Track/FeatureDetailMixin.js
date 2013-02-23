/**
 * Mixin with methods for parsing making default feature detail dialogs.
 */
define("JBrowse/View/Track/FeatureDetailMixin", [
            'dojo/_base/declare',
            'dojo/_base/array',
            'JBrowse/Util',
            'JBrowse/View/FASTA'
        ],
        function( declare, array, Util, FASTAView ) {

return declare(null,{

    _setupEventHandlers: function() {
        // make a default click event handler
        var eventConf = this.config.events || {};
        if( ! eventConf.click ) {
            eventConf.click = (this.config.style||{}).linkTemplate
                    ? { action: "newWindow", url: this.config.style.linkTemplate }
                    : { action: "contentDialog",
                        title: '{type} {name}',
                        content: dojo.hitch( this, 'defaultFeatureDetail' ) };
        }

        // process the configuration to set up our event handlers
        this.eventHandlers = (function() {
            var handlers = dojo.clone( eventConf );
            // find conf vars that set events, like `onClick`
            for( var key in this.config ) {
                var handlerName = key.replace(/^on(?=[A-Z])/, '');
                if( handlerName != key )
                    handlers[ handlerName.toLowerCase() ] = this.config[key];
            }
            // interpret handlers that are just strings to be URLs that should be opened
            for( key in handlers ) {
                if( typeof handlers[key] == 'string' )
                    handlers[key] = { url: handlers[key] };
            }
            return handlers;
        }).call(this);
        this.eventHandlers.click = this._makeClickHandler( this.eventHandlers.click );
    },

    /**
     * Make a default feature detail page for the given feature.
     * @returns {HTMLElement} feature detail page HTML
     */
    defaultFeatureDetail: function( /** JBrowse.Track */ track, /** Object */ f, /** HTMLElement */ featDiv, /** HTMLElement */ container ) {
        var fmt = dojo.hitch( this, '_fmtDetailField' );
        container = container || dojo.create('div', { className: 'detail feature-detail feature-detail-'+track.name, innerHTML: '' } );
        var coreDetails = dojo.create('div', { className: 'core' }, container );
        coreDetails.innerHTML += '<h2 class="sectiontitle">Primary Data</h2>';
        coreDetails.innerHTML += fmt( 'Name', f.get('name') );
        coreDetails.innerHTML += fmt( 'Type', f.get('type') );
        coreDetails.innerHTML += fmt( 'Description', f.get('note') );
        coreDetails.innerHTML += fmt(
            'Position',
            Util.assembleLocString({ start: f.get('start'),
                                     end: f.get('end'),
                                     ref: this.refSeq.name,
                                     strand: f.get('strand')
                                   })
        );
        coreDetails.innerHTML += fmt( 'Length', Util.addCommas(f.get('end')-f.get('start'))+' bp' );

        // render any additional tags as just key/value
        var additionalTags = array.filter( f.tags(), function(t) {
            return ! {name:1,start:1,end:1,strand:1,note:1,subfeatures:1,type:1}[t.toLowerCase()];
        });

        if( additionalTags.length ) {
            var at_html = '<div class="additional"><h2 class="sectiontitle">Attributes</h2>';
            dojo.forEach( additionalTags.sort(), function(t) {
                at_html += fmt( t, f.get(t) );
            });
            at_html += '</div>';
            container.innerHTML += at_html;
        }

        // render the sequence underlying this feature if possible
        var field_container = dojo.create('div', { className: 'field_container feature_sequence' }, container );
        dojo.create( 'h2', { className: 'field feature_sequence', innerHTML: 'Region sequence' }, field_container );
        var valueContainerID = 'feature_sequence'+this._uniqID();
        var valueContainer = dojo.create(
            'div', {
                id: valueContainerID,
                innerHTML: '<div style="height: 12em">Loading...</div>',
                className: 'value feature_sequence'
            }, field_container);
        track.browser.getStore('refseqs', dojo.hitch(this,function( refSeqStore ) {
            valueContainer = dojo.byId(valueContainerID) || valueContainer;
            if( refSeqStore ) {
                refSeqStore.getFeatures(
                    { ref: this.refSeq.name, start: f.get('start'), end: f.get('end')},
                    // feature callback
                    dojo.hitch( this, function( feature ) {
                        var seq = feature.get('seq');
                        valueContainer = dojo.byId(valueContainerID) || valueContainer;
                        valueContainer.innerHTML = '';
                        // the HTML is rewritten by the dojo dialog
                        // parser, but this callback may be called either
                        // before or after that happens.  if the fetch by
                        // ID fails, we have come back before the parse.
                        var textArea = new FASTAView({ width: 62, htmlMaxRows: 10 })
                                           .renderHTML(
                                               { ref:   this.refSeq.name,
                                                 start: f.get('start'),
                                                 end:   f.get('end'),
                                                 strand: f.get('strand'),
                                                 type: f.get('type')
                                               },
                                               f.get('strand') == -1 ? Util.revcom(seq) : seq,
                                               valueContainer
                                           );
                  }),
                  // end callback
                  function() {},
                  // error callback
                  dojo.hitch( this, function() {
                      valueContainer = dojo.byId(valueContainerID) || valueContainer;
                      valueContainer.innerHTML = '<span class="ghosted">reference sequence not available</span>';
                  })
                );
            } else {
                valueContainer.innerHTML = '<span class="ghosted">reference sequence not available</span>';
            }
        }));

        // render any subfeatures this feature has
        var subfeatures = f.get('subfeatures');
        if( subfeatures && subfeatures.length ) {
            this._subfeaturesDetail( track, subfeatures, container );
        }

        return container;
    },

    _uniqID: function() {
        this._idCounter = this._idCounter || 0;
        return this._idCounter++;
    },

    _subfeaturesDetail: function( track, subfeatures, container ) {
            var field_container = dojo.create('div', { className: 'field_container subfeatures' }, container );
            dojo.create( 'h2', { className: 'field subfeatures', innerHTML: 'Subfeatures' }, field_container );
            var subfeaturesContainer = dojo.create( 'div', { className: 'value subfeatures' }, field_container );
            array.forEach( subfeatures || [], function( subfeature ) {
                    this.defaultFeatureDetail(
                        track,
                        subfeature,
                        null,
                        dojo.create('div', {
                                        className: 'detail feature-detail subfeature-detail feature-detail-'+track.name+' subfeature-detail-'+track.name,
                                        innerHTML: ''
                                    }, subfeaturesContainer )
                    );
            },this);
    }

});
});