
// Final Year Project for Flood Mapping using Sentinel Data, completed under the guidance of Dr. G.S Baghel  



// DataSets
var gsw = ee.Image("JRC/GSW1_2/GlobalSurfaceWater"),
    admin2 = ee.FeatureCollection("FAO/GAUL_SIMPLIFIED_500m/2015/level2"),
    hydrosheds = ee.Image("WWF/HydroSHEDS/03VFDEM"),
    imageCollection = ee.ImageCollection("COPERNICUS/S1_GRD");



//Sample Date Range
var beforeStart = '2018-07-15'
var beforeEnd = '2018-08-10'
var afterStart = '2018-08-10'
var afterEnd = '2018-08-23'


// District Data
var ernakulam = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Ernakulam'))
var alappuzha = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Alappuzha'))
var idukki = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Idukki'))
var Kannur = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Kannur'))
var Kasaragod = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Kasaragod'))
var Kollam = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Kollam'))
var Kottayam = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Kottayam'))
var Kozhikode = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Kozhikode'))
var Malappuram = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Malappuram'))
var Palakkad = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Palakkad'))
var Pathanamthitta = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Pathanamthitta'))
var Thiruvananthapuram = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Thiruvananthapuram'))
var Thrissur = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Thrissur'))
var Wayanad = admin2.filter(ee.Filter.eq('ADM2_NAME', 'Wayanad'))


var geometry = ernakulam.geometry()
Map.addLayer(geometry, {color: 'grey'}, ' District')

var collection= ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.eq('instrumentMode','IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')) 
  .filter(ee.Filter.eq('resolution_meters',10))
  .filterBounds(geometry)
  .select('VH');

var beforeCollection = collection.filterDate(beforeStart, beforeEnd)
var afterCollection = collection.filterDate(afterStart,afterEnd)

var before = beforeCollection.mosaic().clip(geometry);
var after = afterCollection.mosaic().clip(geometry);

Map.addLayer(before, {min:-25,max:0}, 'Before Floods', false);
Map.addLayer(after, {min:-25,max:0}, 'After Floods', false); 

var beforeFiltered = ee.Image(toDB(RefinedLee(toNatural(before))))
var afterFiltered = ee.Image(toDB(RefinedLee(toNatural(after))))

Map.addLayer(beforeFiltered, {min:-25,max:0}, 'Before Filtered', false);
Map.addLayer(afterFiltered, {min:-25,max:0}, 'After Filtered', false); 



var histogramBefore = before.select('VH').reduceRegion({
  reducer: ee.Reducer.histogram(255, 2)
      .combine('mean', null, true)
      .combine('variance', null, true), 
  geometry: geometry, 
  scale: 30,
  bestEffort: true
});


print(ui.Chart.image.histogram(before, geometry, 30));



var histogramAfter = after.select('VH').reduceRegion({
  reducer: ee.Reducer.histogram(255, 2)
      .combine('mean', null, true)
      .combine('variance', null, true), 
  geometry: geometry, 
  scale: 30,
  bestEffort: true
});


//print(ui.Chart.image.histogram(after, geometry, 30));


var otsuBeforeThreshold = Otsu(histogramBefore.get('VH_histogram'));
//var otsuAfterThreshold = Otsu(histogramAfter.get('VH_histogram'));

// #############################
// Otsu Threshold
// ############################


function Otsu(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  print(ui.Chart.array.values(ee.Array(bss), 0, means));
  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};


print('Threshold', otsuBeforeThreshold);
var beforeInitial = before.select('VH').lt(otsuBeforeThreshold);
Map.addLayer(beforeInitial.mask(beforeInitial), {palette: 'pink'}, 'Flood OTSU Before', false);


//print('After threshold', otsuAfterThreshold);
// var afterInitial = after.select('VH').lt(otsuAfterThreshold);
// Map.addLayer(afterInitial.mask(afterInitial), {palette: 'pink'}, 'Flood OTSU After', false);


var difference = afterFiltered.divide(beforeFiltered);

// Initial estimate of flooded pixels
var flooded = difference.gt(1.25).rename('water').selfMask();
Map.addLayer(flooded, {min:0, max:1, palette: ['orange']}, 'Initial Flood Area', false);


// Mask out area with permanent/semi-permanent water
var permanentWater = gsw.select('seasonality').gte(5).clip(geometry)
var flooded = flooded.where(permanentWater, 0).selfMask()
Map.addLayer(permanentWater.selfMask(), {min:0, max:1, palette: ['blue']}, 'Permanent Water')

// Mask out areas with more than 5 percent slope using the HydroSHEDS DEM
var slopeThreshold = 5;
var terrain = ee.Algorithms.Terrain(hydrosheds);
var slope = terrain.select('slope').clip(geometry);
var flooded = flooded.updateMask(slope.lt(slopeThreshold));
Map.addLayer(slope.gte(slopeThreshold).selfMask(), {min:0, max:1, palette: ['cyan']}, 'Steep Areas', false)

        
// Remove isolated pixels
// connectedPixelCount is Zoom dependent, so visual result will vary
var connectedPixelThreshold = 8;
var connections = flooded.connectedPixelCount(25)
var flooded = flooded.updateMask(connections.gt(connectedPixelThreshold))
Map.addLayer(connections.lte(connectedPixelThreshold).selfMask(), {min:0, max:1, palette: ['yellow']}, 'Disconnected Areas', false)

Map.addLayer(flooded, {min:0, max:1, palette: ['red']}, 'Flooded Areas');

// Calculate Affected Area
print('Total District Area (Ha)', geometry.area().divide(10000))

var stats = flooded.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: geometry,
  scale: 30,
  maxPixels: 1e10,
  tileScale: 16
})
print('Flooded Area (Ha)', ee.Number(stats.get('water')).divide(10000))


var histogram = flooded.select('VH').reduceRegion({
  reducer: ee.Reducer.histogram(255, 2)
      .combine('mean', null, true)
      .combine('variance', null, true), 
  geometry: geometry, 
  scale: 30,
  bestEffort: true
});


//############################
//
// Speckle Filtering Functions
//
//############################


function RefinedLee(img) {
  // Kernal set up 3x3 
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  //Mean and Variance for Kernal of Matrix
  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use 3x3 windows in 7x7 windows to determine gradient  directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());  

  //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
  //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create necessary stacks for mean and variance using the original kernels.  Mask  with relevant direction.
  var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  //Add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  //Generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
  return(result.arrayFlatten([['sum']]));
}


// Function  utilized to convert to natural
function toNatural(img) {
    return ee.Image(10.0).pow(img.select(0).divide(10.0));
  }
  
  //Function utilized to convert to dB
  function toDB(img) {
    return ee.Image(img).log10().multiply(10.0);
  }
  