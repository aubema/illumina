
 Ext.onReady(function() {
     
     
     //Menu
     var time_cb = new Ext.form.ComboBox({
	    typeAhead: false,
	    id: 'time_cb',
	    width: 220,
	    fieldLabel: 'Simulated time',
	    triggerAction: 'all',
	    lazyRender:true,
	    mode: 'local',
	    store: new Ext.data.ArrayStore({
	        id: 0,
	        fields: [
	            'timeId',
	            'displayText'
	        ],
	        data: [['after_midnight', 'Summer 2009 — After Midnight '], 
	        ['before_midnight', 'Summer 2009 - Before Midnight']]
	    }),
	    value: 'before_midnight',
	    valueField: 'timeId',
	    displayField: 'displayText'
	});
	

var lambda_cb = new Ext.form.ComboBox({
    typeAhead: false,
    id: 'lambda_cb',
    fieldLabel: 'Element - Wavelength',
    triggerAction: 'all',
    lazyRender:true,
    mode: 'local',
    store: new Ext.data.ArrayStore({
        id: 1,
        fields: [
            'lambda',
            'displayText'
        ],
        data: [[435, 'Hg - 435 nm'], [498, 'Na - 498 nm'], [546, 'Hg - 546 nm'],
               [568, 'Na - 568 nm'], [615, 'Na - 615 nm']]
    }),
    value: 615,
    valueField: 'lambda',
    displayField: 'displayText'
});

var site_cb = new Ext.form.ComboBox({
	typeAhead: false,
	width: 220,
	id: 'site_cb',
	fieldLabel: 'Site',
	triggerAction: 'all',
	mode: 'local',
	store: new Ext.data.ArrayStore({
		id: 2,
		fields: [
			'site_abr',
			'displayText'
		],
		data: [
			['OT', 'Observatorio del Teide'], 
			['ORM', 'Observatorio Roque de los Muchachos']
		]
	}),
	value: 'OT',
	valueField: 'site_abr',
	displayField: 'displayText'
});

var site_section = new Ext.form.FieldSet({
    title: 'Simulation selection',
    labelWidth: 90, // label settings here cascade unless overridden
    autoHeight: true,
    layout: 'form',
    frame:true,
    bodyStyle:'padding:5px 5px 0',
    items: [site_cb, time_cb]
});

var map_type_section = new Ext.form.FieldSet({
	title: 'Map type',
    labelWidth: 5, // label settings here cascade unless overridden
    autoHeight: true,
    layout: 'form',
    frame:true,
    bodyStyle:'padding:5px 5px 0',
    defaults: {      // defaults applied to items
        layout: 'form',
        border: false,
        bodyStyle: 'padding:4px'
    },
    items: {xtype: 'radiogroup',
        	fieldLabel: '',
        	id: 'type_rd',
        	width: 400,
        	autoHeight: true,
        	items:
	        [{
	      		type: 'radio',
	      		inputValue: 'PCL',
	      		checked: true,
	      		name: 'map_type',
	      		boxLabel: 'Relative Contribution Map'
	        },{
	        	type: 'radio',
	      		inputValue: 'PCW',
	      		name: 'map_type',
	      		boxLabel: 'Relative Sensivity Map'
	        }]
	        }
});


var aod_cb = new Ext.form.ComboBox({
	typeAhead: false,
	id: 'aod_cb',
	labelWidth: 160,
	triggerAction: 'all',
	fieldLabel: 'Aerosol optical depth',
	mode: 'local',
	store: new Ext.data.ArrayStore({
		id: 3,
		fields: [
			'aod_val',
			'displayText'
		],
		data: [
			['ta0.025', '0.025'], 
			['ta0.050', '0.050'],
			['ta0.100', '0.100'], 
			['ta0.200', '0.200'],
		]
	}),
	value: 'ta0.200',
	valueField: 'aod_val',
	displayField: 'displayText'
});


var angles_el = new Ext.form.ComboBox({
	typeAhead: false,
	id: 'angles_el',
	triggerAction: 'all',
	fieldLabel: 'Elevation angle (deg.)',
	mode: 'local',
	store: new Ext.data.ArrayStore({
		id: 4,
		fields: ['anglesel', 'displayText'],
		data: [
		    ['el5', '5'],
	            ['el10', '10'],
                    ['el20', '20'],
                    ['el30', '30'],
                    ['el45', '45'],
                    ['el70', '70'],
                    ['el90', '90'],
		]
	}),
	value: 'el90',
	valueField: 'anglesel',
	displayField: 'displayText'
});


var angles_az = new Ext.form.ComboBox({
        typeAhead: false,
        id: 'angles_az',
        triggerAction: 'all',
        fieldLabel: 'Azimut angle (deg.) N=0, E=90...',
        mode: 'local',
        store: new Ext.data.ArrayStore({
                id: 4,
                fields: ['anglesaz', 'displayText'],
                data: [
                    ['az0', '0'],
                    ['az22', '22'],
                    ['az30', '30'],
                    ['az45', '45'],
                    ['az60', '60'],
                    ['az67', '67'],
                    ['az90', '90'],
                    ['az112', '112'],
                    ['az120', '120'],
                    ['az135', '135'],
                    ['az150', '150'],
                    ['az157', '157'],
                    ['az180', '180'],
                    ['az202', '202'],
                    ['az210', '210'],
                    ['az225', '225'],
                    ['az240', '240'],
                    ['az257', '257'],
                    ['az270', '270'],
                    ['az292', '292'],
                    ['az300', '300'],
                    ['az315', '315'],
                    ['az330', '330'],
                    ['az337', '337'],
                ]
        }),
        value: 'az0',
        valueField: 'anglesaz',
        displayField: 'displayText'
});

var physical_section = new Ext.form.FieldSet({
	title: 'Physical parameters',
	collapsed: false,
	collapsible: true,
    labelWidth: 160, // label settings here cascade unless overridden
    autoHeight: true,
    layout: 'form',
    frame:true,
    bodyStyle:'padding:5px 5px 0',
    defaults: {      // defaults applied to items
        layout: 'form',
        border: false,
        bodyStyle: 'padding:4px'
    },
    items: [lambda_cb, 
          aod_cb,
          angles_el,
          angles_az]
       });


    var simple = new Ext.FormPanel({
        id: 'param_form',
        labelWidth: 180,
        frame: true,
        title: 'Parameter',
        bodyStyle:'padding:5px 5px 0',
        width: 800,

        items: [site_section,
          map_type_section,
          physical_section
        ],

        buttons: [{
            text: 'Submit',
        handler: function() {
         var c_map_type = Ext.getCmp('param_form').findById('type_rd').getValue().inputValue;
         var c_site = Ext.getCmp('param_form').findById('site_cb').getValue(asString=true);
         var c_aod = Ext.getCmp('param_form').findById('aod_cb').getValue(asString=true);
         var c_time = Ext.getCmp('param_form').findById('time_cb').getValue(asString=true);
         var c_lambda = Ext.getCmp('param_form').findById('lambda_cb').getValue(asString=true);
         var c_anglesel =  Ext.getCmp('param_form').findById('angles_el').getValue(asString=true);
         var c_anglesaz =  Ext.getCmp('param_form').findById('angles_az').getValue(asString=true);
         if (angles_el = "el90")
         {
         c_anglesaz="az0"
         }
         var layerName = c_map_type + '-' + c_site + '-' + c_time + '-rd4000-' + c_aod + '-wl' + c_lambda + '-' + c_anglesel + '-' + c_anglesaz + '_epsg900913';
         mappanel = Ext.getCmp('main_mappanel');
         layerO = mappanel.map.getLayersByName("Light Pollution Atlas Current Map")[0];
         mappanel.map.removeLayer(layerO);
         var lp_layer = new OpenLayers.Layer.WMS( "Light Pollution Atlas Current Map",
                    "http://galileo.graphycs.cegepsherbrooke.qc.ca/lpawms2", 
                    {layers: layerName, 
                    transparent:true},
                    {isBaseLayer: false, visibility: true, alpha: true}  );
         mappanel.map.addLayer(lp_layer);
        }
        },{
            text: 'Reset'
        }],
        applyTo: 'menu'
    });

physical_section.collapse()
     // Map
     
        var map 				= new OpenLayers.Map("map", {projection:'EPSG:900913'});
        arrayOSM = ["http://otile1.mqcdn.com/tiles/1.0.0/osm/${z}/${x}/${y}.jpg",
                    "http://otile2.mqcdn.com/tiles/1.0.0/osm/${z}/${x}/${y}.jpg",
                    "http://otile3.mqcdn.com/tiles/1.0.0/osm/${z}/${x}/${y}.jpg",
                    "http://otile4.mqcdn.com/tiles/1.0.0/osm/${z}/${x}/${y}.jpg"];
        arrayAerial = ["http://oatile1.mqcdn.com/tiles/1.0.0/sat/${z}/${x}/${y}.jpg",   
                        "http://oatile2.mqcdn.com/tiles/1.0.0/sat/${z}/${x}/${y}.jpg",
                        "http://oatile3.mqcdn.com/tiles/1.0.0/sat/${z}/${x}/${y}.jpg",
                        "http://oatile4.mqcdn.com/tiles/1.0.0/sat/${z}/${x}/${y}.jpg"];
        var osm 			= new OpenLayers.Layer.OSM("MaQuest-OSM Tiles", arrayOSM);
        var osm_sat         = new OpenLayers.Layer.OSM("MapQuest Open Aerial Tiles", arrayAerial);
        var google_l        = new OpenLayers.Layer.Google("Google Sat", {type: google.maps.MapTypeId.SATELLITE, numZoomLevels: 22});
        var fromProjection = new OpenLayers.Projection("EPSG:4326");   // Transform from WGS 1984
        var toProjection   = new OpenLayers.Projection("EPSG:900913"); // to Spherical Mercator Projection
        var position       = new OpenLayers.LonLat(-16.62,28.28).transform( fromProjection, toProjection);
        var init_lp_layer = 'PCL-OT-before_midnight-rd4000-ta0.200-wl615-el5-az60_epsg900913';
        var lp_layer = new OpenLayers.Layer.WMS( "Light Pollution Atlas Current Map",
                    "http://galileo.graphycs.cegepsherbrooke.qc.ca/lpawms2", 
                    {layers: init_lp_layer, 
                    transparent:true},
                    {isBaseLayer: false, visibility: true, alpha: true}  );
        
        map.addLayers([osm, osm_sat, google_l, lp_layer]);
        map.addControl(new OpenLayers.Control.LayerSwitcher());

         mapPanel = new GeoExt.MapPanel({
         	id: 'main_mappanel',
            height: 400,
            width: 600,
            region: 'center',
            map: map,
            title: 'Light Pollution contribution or sensitivity Map',
            center: position,
            zoom: 6
        });
        
        legendPanel = new GeoExt.LegendPanel({
	        defaults: {
	            labelCls: 'mylabel',
	            style: 'padding:5px'
	        },
	        bodyStyle: 'padding:5px',
	        width: 220,
	        autoScroll: true,
	        region: 'west'




    	});
    	
    	new Ext.Panel({
	        title: "Map",
	        layout: 'border',
	        renderTo: 'gxmap',
	        height: 500,
	        width: 800,
	        items: [legendPanel, mapPanel]
	    });
    });
