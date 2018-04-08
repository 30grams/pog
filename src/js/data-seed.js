import SimplexNoise from './vendor/simplex-noise.js';

'use strict';

const mapWidth = 20,
  mapHeight = 30;

// HASHING 
String.prototype.hashCode = function() {
  // Java like string Hash function 
  var hash = 0, i, chr, len;
  if (this.length == 0) return hash;
  for (i = 0, len = this.length; i < len; i++) {
    chr   = this.charCodeAt(i);
    hash  = ((hash << 5) - hash) + chr;
    hash |= 0; // Convert to 32bit integer
  }
  return hash;
};

// SEEDING
export function generateMap(seed) {
  if (!seed) {
    seed = Math.random.toString();
  }

  const simplex = new SimplexNoise(seed.hashCode());
  // initializing a new simplex instance
  // do this only once as it is relatively expensive

  let assignElevation = (x,y) => {
    return simplex.noise2D(x,y); // just for now
  }

  let elevationMap = new Array(mapWidth).fill(null).map((el, i) => {
    return new Array(mapHeight).fill(null).map((el, j) => {
      // run each object through the constructor
      return assignElevation(i, j);
    })
  });

  return elevationMap;
}


// BIOMES
function biome(altitude, temperature, moisture) {
  return 'tbd';
}


