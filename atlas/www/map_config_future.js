
 Ext.onReady(function() {
     
     
     //Menu

var sitetime_cb = new Ext.form.ComboBox({
        typeAhead: false,
        width: 420,
        id: 'sitetime_cb',
        fieldLabel: 'Site & time',
        triggerAction: 'all',
        mode: 'local',
        store: new Ext.data.ArrayStore({
                id: 2,
                fields: [
                        'sitetime_abr',
                        'displayText'
                ],
                data: [
                        ['OT-b_midnight_winter_2009', 'Observatorio del Teide - before midnight, winter 2009' ],
                        ['OT-a_midnight_winter_2009', 'Observatorio del Teide - after midnight, winter 2009' ],
                        ['ORM-b_midnight_winter_2009', 'Observatorio Roque de los Muchachos - before midnight, winter 2009'],
                        ['ORM-a_midnight_winter_2009', 'Observatorio Roque de los Muchachos - after midnight, winter 2009'],
                        ['OMM-b_midnight_summer_2005', 'Observatoire du Mont-Megantic - before midnight, summer 2005'],
                        ['OMM-b_midnight_winter_2005', 'Observatoire du Mont-Megantic - before midnight, winter 2005'],
                        ['OMM-b_midnight_summer_2009', 'Observatoire du Mont-Megantic - before midnight, summer 2009'],
                        ['OMM-b_midnight_winter_2009', 'Observatoire du Mont-Megantic - before midnight, winter 2009']
                ]
        }),
        value: 'OT-b_midnight_winter_2009',
        valueField: 'sitetime_abr',
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
        data: [[436, 'Hg - 436 nm'], [498, 'Na - 498 nm'], [546, 'Hg - 546 nm'],
               [569, 'Na - 569 nm'], [616, 'Na - 616 nm']]
    }),
    value: 569,
    valueField: 'lambda',
    displayField: 'displayText'
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
			['ta0.05', '0.05'], 
			['ta0.1', '0.1'],
			['ta0.2', '0.2'], 
			['ta0.5', '0.5'],
                        ['ta1.0', '1.0'],
		]
	}),
	value: 'ta0.1',
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
		    ['el15', '15'],
                    ['el20', '20'],
                    ['el30', '30'],
                    ['el40', '40'],
                    ['el60', '60'],
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
                id: 5,
                fields: ['anglesaz', 'displayText'],
                data: [
                    ['az0', '0'],
                    ['az15', '15'],
                    ['az30', '30'],
                    ['az45', '45'],
                    ['az60', '60'],
                    ['az75', '75'],
                    ['az90', '90'],
                    ['az105', '105'],
                    ['az120', '120'],
                    ['az135', '135'],
                    ['az150', '150'],
                    ['az165', '165'],
                    ['az180', '180'],
                    ['az195', '195'],
                    ['az210', '210'],
                    ['az225', '225'],
                    ['az240', '240'],
                    ['az255', '255'],
                    ['az270', '270'],
                    ['az285', '285'],
                    ['az300', '300'],
                    ['az315', '315'],
                    ['az330', '330'],
                    ['az345', '345'],
                ]
        }),
        value: 'az0',
        valueField: 'anglesaz',
        displayField: 'displayText'
});



    var simple = new Ext.FormPanel({
        id: 'param_form',
        labelWidth: 180,
        frame: true,
        title: 'Parameter',
        bodyStyle:'padding:5px 5px 0',
        width: 800,

        items: [{xtype: 'radiogroup',
        	fieldLabel: 'Map type',
        	id: 'type_rd',
        	width: 400,
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
	        },
          sitetime_cb,
          lambda_cb, 
          aod_cb,
          angles_el,
          angles_az
        ],

        buttons: [{
            text: 'Submit',
        handler: function() {
         var c_map_type = Ext.getCmp('param_form').findById('type_rd').getValue().inputValue;
         var c_sitetime = Ext.getCmp('param_form').findById('sitetime_cb').getValue(asString=true);
         var c_aod = Ext.getCmp('param_form').findById('aod_cb').getValue(asString=true);
         var c_lambda = Ext.getCmp('param_form').findById('lambda_cb').getValue(asString=true);
         var c_anglesel =  Ext.getCmp('param_form').findById('angles_el').getValue(asString=true);
         var c_anglesaz =  Ext.getCmp('param_form').findById('angles_az').getValue(asString=true);
         if (angles_el == "el90")
         {
         c_anglesaz="az0"
         }
         var layerName = c_map_type + '-' + c_sitetime + '-' + c_aod + '-wl' + c_lambda + '-' + c_anglesel + '-' + c_anglesaz + '_epsg900913';
         mappanel = Ext.getCmp('main_mappanel');
         layerO = mappanel.map.getLayersByName("Light Pollution Atlas Current Map")[0];
         mappanel.map.removeLayer(layerO);
         var lp_layer = new OpenLayers.Layer.WMS( "Light Pollution Atlas Current Map",
                    "http://galileo.graphycs.cegepsherbrooke.qc.ca/lpawms3", 
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
        var init_lp_layer = 'PCL-OT-b_midnight_2009-ta0.200-wl615-el5-az60_epsg900913';
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
