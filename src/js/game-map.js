import { CylinderGeometry , Mesh, MeshLambertMaterial } from 'three';

const tileRadius = 1, 
 tileDepth = 5;
 // represents cyliner size in meters

const tileHeight = tileRadius * 2;
const tileWidth = Math.sqrt(3)/2 * tileHeight;

export default function populate(scene) {

  // TILES
  let tileGeometry = new CylinderGeometry( tileRadius, tileRadius, tileDepth, 6, 1 );
  let material = new MeshLambertMaterial( { color: 0x77dd77, flatShading: true, overdraw: 0.5} );

  tileGeometry.computeFlatVertexNormals();


  for ( var i = -2; i < 3; i += 1 ) {
    for ( var j = -2; j < 3; j += 1 ) {
      // axial coordinates

      var tile = new Mesh( tileGeometry, material );
      tile.position.x = tileWidth * ( i + j / 2 ); 
      tile.position.z = tileHeight * 3/4 * j; 
      tile.position.y = ( -Math.floor(Math.random() * 5) * 0.2) - tileDepth/2; 

      tile.castShadow = true;
      tile.receiveShadow = true;
      scene.add( tile );
    }
  }
}
