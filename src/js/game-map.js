import { CylinderGeometry , Mesh, MeshLambertMaterial } from 'three';

const tileWidth = 1,
 tileHeight = 5;
const tileRadius = tileWidth / Math.sqrt(3);

export default function populate(scene) {

  // TILES
  let tileGeometry = new CylinderGeometry( tileRadius, tileRadius, tileHeight, 6, 1 );
  let material = new MeshLambertMaterial( { color: 0x77dd77, flatShading: true, overdraw: 0.5} );

  tileGeometry.computeFlatVertexNormals();


  for ( var i = -2; i < 3; i += 1 ) {
    for ( var j = -2; j < 3; j += 1 ) {

      var tile = new Mesh( tileGeometry, material );
      tile.position.x = tileWidth * ( i + j * 0.5 ); 
      tile.position.z = (tileWidth * 0.5 *  Math.sqrt(3)) * j; 
      tile.position.y = ( -Math.floor(Math.random() * 5) * 0.2) - tileHeight/2; 

      tile.castShadow = true;
      tile.receiveShadow = true;
      scene.add( tile );
    }
  }
}
