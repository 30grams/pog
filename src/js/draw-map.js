import { CylinderGeometry , Mesh, MeshLambertMaterial } from 'three';
import {generateMap} from './data-seed.js';
import {flatMaterial} from './draw-textures.js';

const tileRadius = 1, 
 tileDepth = 5;
 // represents cyliner size in meters

const tileHeight = tileRadius * 2;
const tileWidth = Math.sqrt(3)/2 * tileHeight;

export default function populate(scene) {

  // TILES
  let tileGeometry = new CylinderGeometry( tileRadius, tileRadius, tileDepth, 6, 1 );
  let material = new MeshLambertMaterial( { color: 0x77dd77, flatShading: true} );

  tileGeometry.computeFlatVertexNormals();

  // MAP
  let tileMap = generateMap('pepe');

/*
  tileMap.forEach(function(el, i) {
    el.forEach(function(el, j) {
      // axial coordinates

      var tile = new Mesh( tileGeometry, material );
      tile.position.x = tileWidth * ( i + j / 2 ); 
      tile.position.z = tileHeight * 3/4 * j; 
      tile.position.y = ( -Math.floor(Math.random() * 5) * 0.2) - tileDepth/2; 

      tile.castShadow = true;
      tile.receiveShadow = true;
      scene.add( tile );
    })
  });
*/

  for ( var i = -Math.round(tileMap.length / 2); i < tileMap.length - Math.round(tileMap.length / 2); i++ ) {
    for ( var j = -Math.round(tileMap[0].length / 2); j < tileMap[0].length - Math.round(tileMap[0].length / 2); j++ ) {
      // axial coordinates

      var tile = new Mesh( tileGeometry, flatMaterial() );
      tile.position.x = tileWidth * ( i + j / 2 ); 
      tile.position.z = tileHeight * 3/4 * j; 
      tile.position.y = ( -Math.floor(Math.random() * 5) * 0.2) - tileDepth/2; 

      tile.castShadow = true;
      tile.receiveShadow = true;
      scene.add( tile );
    }
  }
  
}
