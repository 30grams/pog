/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */
'use strict';

const hashInt = require('hash-int');
/*

'use strict';

var A
if(typeof Uint32Array === undefined) {
  A = [ 0 ]
} else {
  A = new Uint32Array(1)
}

function hashInt(x) {
  A[0]  = x|0
  A[0] -= (A[0]<<6)
  A[0] ^= (A[0]>>>17)
  A[0] -= (A[0]<<9)
  A[0] ^= (A[0]<<4)
  A[0] -= (A[0]<<3)
  A[0] ^= (A[0]<<10)
  A[0] ^= (A[0]>>>15)
  return A[0]


*/


/*


String.prototype.hashCode = function() {
  var hash = 0, i, chr;
  if (this.length === 0) return hash;
  for (i = 0; i < this.length; i++) {
    chr   = this.charCodeAt(i);
    hash  = ((hash << 5) - hash) + chr;
    hash |= 0; // Convert to 32bit integer
  }
  return hash;
};

//https://stackoverflow.com/a/7616484
*/


/*
var mystring = "t7";

function pad (str, max) {
  str = str.toString();
  return str.length < max ? pad("0" + str, max) : str;
}


String.prototype.hashCode = function() {
  var hash = 0, i, chr, len;
  if (this.length == 0) return hash;
  for (i = 0, len = this.length; i < len; i++) {
    chr   = this.charCodeAt(i);
    hash  = ((hash << 5) - hash) + chr;
    hash |= 0; // Convert to 32bit integer
  }
  return hash;
};

var max_length =7;
$('#result').text(mystring.hashCode());
$('#resultpad').text(pad(mystring.hashCode(), max_length));     //just padded
$('#resultpadsubstr').text(pad(mystring.hashCode().toString().substr(0,max_length), max_length)); //padded and cut above 7

// https://stackoverflow.com/a/33607807
*/


exports.makeRandInt = function(seed) {
    let i = 0;
    return function(N) {
        i++;
        return hashInt(seed + i) % N;
    };
};

exports.makeRandFloat = function(seed) {
    let randInt = exports.makeRandInt(seed);
    let divisor = 0x10000000;
    return function() {
        return randInt(divisor) / divisor;
    };
};







