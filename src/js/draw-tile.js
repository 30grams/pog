import Mesh from 'three';


export tileMesh = function() {
    // Run the Mesh constructor with the given arguments
    Mesh.apply(this, arguments);
};

// Make tileMesh have the same methods as Mesh
tileMesh.prototype = Object.create(Mesh.prototype);
// Make sure the right constructor gets called
tileMesh.prototype.constructor = tileMesh;

// HOVER
tileMesh.prototype.mouseover = function(){
  this.material.color.setHex( 0xffff00 )
};



