(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
/*
 * From https://github.com/redblobgames/dual-mesh
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 *
 * Generate a random triangle mesh for the area 0 <= x <= 1000, 0 <= y <= 1000
 *
 * This program runs on the command line (node)
 */

'use strict';

let Poisson = require('poisson-disk-sampling'); // MIT licensed
let Delaunator = require('delaunator');        // ISC licensed

function s_next_s(s) { return (s % 3 == 2) ? s-2 : s+1; }


function checkPointInequality({_r_vertex, _s_start_r, _s_opposite_s}) {
    // TODO: check for collinear vertices. Around each red point P if
    // there's a point Q and R both connected to it, and the angle P→Q and
    // the angle P→R are 180° apart, then there's collinearity. This would
    // indicate an issue with poisson disc point selection.
}


function checkTriangleInequality({_r_vertex, _s_start_r, _s_opposite_s}) {
    // check for skinny triangles
    const badAngleLimit = 30;
    let summary = new Array(badAngleLimit).fill(0);
    let count = 0;
    for (let s = 0; s < _s_start_r.length; s++) {
        let r0 = _s_start_r[s],
            r1 = _s_start_r[s_next_s(s)],
            r2 = _s_start_r[s_next_s(s_next_s(s))];
        let p0 = _r_vertex[r0],
            p1 = _r_vertex[r1],
            p2 = _r_vertex[r2];
        let d0 = [p0[0]-p1[0], p0[1]-p1[1]];
        let d2 = [p2[0]-p1[0], p2[1]-p1[1]];
        let dotProduct = d0[0] * d2[0] + d0[1] + d2[1];
        let angleDegrees = 180 / Math.PI * Math.acos(dotProduct);
        if (angleDegrees < badAngleLimit) {
            summary[angleDegrees|0]++;
            count++;
        }
    }
    // NOTE: a much faster test would be the ratio of the inradius to
    // the circumradius, but as I'm generating these offline, I'm not
    // worried about speed right now
    
    // TODO: consider adding circumcenters of skinny triangles to the point set
    if (createMesh.PRINT_WARNINGS && count > 0) {
        console.log('  bad angles:', summary.join(" "));
    }
}


function checkMeshConnectivity({_r_vertex, _s_start_r, _s_opposite_s}) {
    // 1. make sure each side's opposite is back to itself
    // 2. make sure region-circulating starting from each side works
    let ghost_r = _r_vertex.length - 1, out_s = [];
    for (let s0 = 0; s0 < _s_start_r.length; s0++) {
        if (_s_opposite_s[_s_opposite_s[s0]] !== s0) {
            console.log(`FAIL _s_opposite_s[_s_opposite_s[${s0}]] !== ${s0}`);
        }
        let s = s0, count = 0;
        out_s.length = 0;
        do {
            count++; out_s.push(s);
            s = s_next_s(_s_opposite_s[s]);
            if (count > 100 && _s_start_r[s0] !== ghost_r) {
                console.log(`FAIL to circulate around region with start side=${s0} from region ${_s_start_r[s0]} to ${_s_start_r[s_next_s(s0)]}, out_s=${out_s}`);
                break;
            }
        } while (s !== s0);
    }
}


/*
 * Add vertices evenly along the boundary of the mesh;
 * use a slight curve so that the Delaunay triangulation
 * doesn't make long thing triangles along the boundary.
 * These points also prevent the Poisson disc generator
 * from making uneven points near the boundary.
 */
function addBoundaryPoints(spacing, size) {
    let N = Math.ceil(size/spacing);
    let points = [];
    for (let i = 0; i <= N; i++) {
        let t = (i + 0.5) / (N + 1);
        let w = size * t;
        let offset = Math.pow(t - 0.5, 2);
        points.push([offset, w], [size-offset, w]);
        points.push([w, offset], [w, size-offset]);
    }
    return points;
}


function addGhostStructure({_r_vertex, _s_start_r, _s_opposite_s}) {
    const numSolidSides = _s_start_r.length;
    const ghost_r = _r_vertex.length;
    
    let numUnpairedSides = 0, firstUnpairedEdge = -1;
    let r_unpaired_s = []; // seed to side
    for (let s = 0; s < numSolidSides; s++) {
        if (_s_opposite_s[s] === -1) {
            numUnpairedSides++;
            r_unpaired_s[_s_start_r[s]] = s;
            firstUnpairedEdge = s;
        }
    }

    let r_newvertex = _r_vertex.concat([[500, 500]]);
    let s_newstart_r = new Int32Array(numSolidSides + 3 * numUnpairedSides);
    s_newstart_r.set(_s_start_r);
    let s_newopposite_s = new Int32Array(numSolidSides + 3 * numUnpairedSides);
    s_newopposite_s.set(_s_opposite_s);

    for (let i = 0, s = firstUnpairedEdge;
         i < numUnpairedSides;
         i++, s = r_unpaired_s[s_newstart_r[s_next_s(s)]]) {

        // Construct a ghost side for s
        let ghost_s = numSolidSides + 3 * i;
        s_newopposite_s[s] = ghost_s;
        s_newopposite_s[ghost_s] = s;
        s_newstart_r[ghost_s] = s_newstart_r[s_next_s(s)];
        
        // Construct the rest of the ghost triangle
        s_newstart_r[ghost_s + 1] = s_newstart_r[s];
        s_newstart_r[ghost_s + 2] = ghost_r;
        let k = numSolidSides + (3 * i + 4) % (3 * numUnpairedSides);
        s_newopposite_s[ghost_s + 2] = k;
        s_newopposite_s[k] = ghost_s + 2;
    }

    return {
        numSolidSides,
        _r_vertex: r_newvertex,
        _s_start_r: s_newstart_r,
        _s_opposite_s: s_newopposite_s
    };
}


/**
 * Create mesh data in a 1000x1000 space for passing to DualMesh
 *
 * Either pass {spacing, random} to choose random spaced points with
 * boundary points, or pass {points} to use an existing set of points
 * without added boundary points. Or pass {spacing, points, random} to 
 * use a given set of points, plus random spaced points, plus boundary
 * points.
 *
 * The mesh generator runs some sanity checks but does not correct the
 * generated points.
 *
 * This interface is insufficient to cover all the possible variants
 * so it is SUBJECT TO CHANGE.
 */
function createMesh({spacing=Infinity, points=[], random=Math.random}) {
    let generator = new Poisson([1000, 1000], spacing, undefined, undefined, random);
    let boundaryPoints = isFinite(spacing)? addBoundaryPoints(spacing, 1000) : [];
    boundaryPoints.forEach((p) => generator.addPoint(p));
    points.forEach((p) => generator.addPoint(p));
    let allPoints = generator.fill();

    let delaunator = new Delaunator(allPoints);
    let graph = {
        _r_vertex: allPoints,
        _s_start_r: delaunator.triangles,
        _s_opposite_s: delaunator.halfedges
    };

    checkPointInequality(graph);
    checkTriangleInequality(graph);
    
    graph = addGhostStructure(graph);
    graph.numBoundaryRegions = boundaryPoints.length;
    checkMeshConnectivity(graph);

    return graph;
}


createMesh.PRINT_WARNINGS = false;
module.exports = createMesh;

},{"delaunator":3,"poisson-disk-sampling":8}],2:[function(require,module,exports){
/*
 * From https://github.com/redblobgames/dual-mesh
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */

'use strict';

/**
 * Represent a triangle-polygon dual mesh with:
 *   - Regions (r)
 *   - Sides (s)
 *   - Triangles (t)
 *
 * Each element has an id:
 *   - 0 <= r < numRegions
 *   - 0 <= s < numSides
 *   - 0 <= t < numTriangles
 *
 * Naming convention: x_name_y takes x (r, s, t) as input and produces
 * y (r, s, t) as output. If the output isn't a mesh index (r, s, t)
 * then the _y suffix is omitted.
 *
 * A side is directed. If two triangles t0, t1 are adjacent, there will
 * be two sides representing the boundary, one for t0 and one for t1. These
 * can be accessed with s_inner_t and s_outer_t.
 *
 * A side also represents the boundary between two regions. If two regions
 * r0, r1 are adjacent, there will be two sides representing the boundary,
 * s_begin_r and s_end_r.
 *
 * Each side will have a pair, accessed with s_opposite_s.
 *
 * The mesh has no boundaries; it wraps around the "back" using a
 * "ghost" region. Some regions are marked as the boundary; these are
 * connected to the ghost region. Ghost triangles and ghost sides
 * connect these boundary regions to the ghost region. Elements that
 * aren't "ghost" are called "solid".
 */
class TriangleMesh {
    static s_to_t(s)   { return (s/3) | 0; }
    static s_prev_s(s) { return (s % 3 == 0) ? s+2 : s-1; }
    static s_next_s(s) { return (s % 3 == 2) ? s-2 : s+1; }

    /**
     * constructor takes partial mesh information and fills in the rest; the
     * partial information is generated in create.js or in deserialize.js
     */
    constructor ({numBoundaryRegions, numSolidSides, _r_vertex, _s_start_r, _s_opposite_s}) {
        Object.assign(this, {numBoundaryRegions, numSolidSides,
                             _r_vertex, _s_start_r, _s_opposite_s});

        this.numSides = _s_start_r.length;
        this.numRegions = _r_vertex.length;
        this.numSolidRegions = this.numRegions - 1;
        this.numTriangles = this.numSides / 3;
        this.numSolidTriangles = this.numSolidSides / 3;
        
        // Construct an index for finding sides connected to a region
        this._r_any_s = new Int32Array(this.numRegions);
        for (let s = 0; s < _s_start_r.length; s++) {
            this._r_any_s[_s_start_r[s]] = this._r_any_s[_s_start_r[s]] || s;
        }

        // Construct triangle coordinates
        this._t_vertex = new Array(this.numTriangles);
        for (let s = 0; s < _s_start_r.length; s += 3) {
            let a = _r_vertex[_s_start_r[s]],
                b = _r_vertex[_s_start_r[s+1]],
                c = _r_vertex[_s_start_r[s+2]];
            if (this.s_ghost(s)) {
                // ghost triangle center is just outside the unpaired side
                let dx = b[0]-a[0], dy = b[1]-a[1];
                this._t_vertex[s/3] = [a[0] + 0.5*(dx+dy), a[1] + 0.5*(dy-dx)];
            } else {
                // solid triangle center is at the centroid
                this._t_vertex[s/3] = [(a[0] + b[0] + c[0])/3,
                                     (a[1] + b[1] + c[1])/3];
            }
        }
    }

    r_x(r)        { return this._r_vertex[r][0]; }
    r_y(r)        { return this._r_vertex[r][1]; }
    t_x(r)        { return this._t_vertex[r][0]; }
    t_y(r)        { return this._t_vertex[r][1]; }
    r_pos(out, r) { out.length = 2; out[0] = this.r_x(r); out[1] = this.r_y(r); return out; }
    t_pos(out, t) { out.length = 2; out[0] = this.t_x(t); out[1] = this.t_y(t); return out; }
    
    s_begin_r(s)  { return this._s_start_r[s]; }
    s_end_r(s)    { return this._s_start_r[TriangleMesh.s_next_s(s)]; }

    s_inner_t(s)  { return TriangleMesh.s_to_t(s); }
    s_outer_t(s)  { return TriangleMesh.s_to_t(this._s_opposite_s[s]); }

    s_next_s(s)   { return TriangleMesh.s_next_s(s); }
    s_prev_s(s)   { return TriangleMesh.s_prev_s(s); }
    
    s_opposite_s(s) { return this._s_opposite_s[s]; }
    
    t_circulate_s(out_s, t) { out_s.length = 3; for (let i = 0; i < 3; i++) { out_s[i] = 3*t + i; } return out_s; }
    t_circulate_r(out_r, t) { out_r.length = 3; for (let i = 0; i < 3; i++) { out_r[i] = this._s_start_r[3*t+i]; } return out_r; }
    t_circulate_t(out_t, t) { out_t.length = 3; for (let i = 0; i < 3; i++) { out_t[i] = this.s_outer_t(3*t+i); } return out_t; }
    
    r_circulate_s(out_s, r) {
        const s0 = this._r_any_s[r];
        let s = s0;
        out_s.length = 0;
        do {
            out_s.push(s);
            s = TriangleMesh.s_next_s(this._s_opposite_s[s]);
        } while (s != s0);
        return out_s;
    }

    r_circulate_r(out_r, r) {
        const s0 = this._r_any_s[r];
        let s = s0;
        out_r.length = 0;
        do {
            out_r.push(this.s_end_r(s));
            s = TriangleMesh.s_next_s(this._s_opposite_s[s]);
        } while (s != s0);
        return out_r;
    }
    
    r_circulate_t(out_t, r) {
        const s0 = this._r_any_s[r];
        let s = s0;
        out_t.length = 0;
        do {
            out_t.push(TriangleMesh.s_to_t(s));
            s = TriangleMesh.s_next_s(this._s_opposite_s[s]);
        } while (s != s0);
        return out_t;
    }

    ghost_r()     { return this.numRegions - 1; }
    s_ghost(s)    { return s >= this.numSolidSides; }
    r_ghost(r)    { return r == this.numRegions - 1; }
    t_ghost(t)    { return this.s_ghost(3 * t); }
    s_boundary(s) { return this.s_ghost(s) && (s % 3 == 0); }
    r_boundary(r) { return r < this.numBoundaryRegions; }
}

module.exports = TriangleMesh;

},{}],3:[function(require,module,exports){
'use strict';

module.exports = Delaunator;

function Delaunator(points, getX, getY) {

    if (!getX) getX = defaultGetX;
    if (!getY) getY = defaultGetY;

    var minX = Infinity;
    var minY = Infinity;
    var maxX = -Infinity;
    var maxY = -Infinity;

    var coords = this.coords = [];
    var ids = this.ids = new Uint32Array(points.length);

    for (var i = 0; i < points.length; i++) {
        var p = points[i];
        var x = getX(p);
        var y = getY(p);
        ids[i] = i;
        coords[2 * i] = x;
        coords[2 * i + 1] = y;
        if (x < minX) minX = x;
        if (y < minY) minY = y;
        if (x > maxX) maxX = x;
        if (y > maxY) maxY = y;
    }

    var cx = (minX + maxX) / 2;
    var cy = (minY + maxY) / 2;

    var minDist = Infinity;
    var i0, i1, i2;

    // pick a seed point close to the centroid
    for (i = 0; i < points.length; i++) {
        var d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
        if (d < minDist) {
            i0 = i;
            minDist = d;
        }
    }

    minDist = Infinity;

    // find the point closest to the seed
    for (i = 0; i < points.length; i++) {
        if (i === i0) continue;
        d = dist(coords[2 * i0], coords[2 * i0 + 1], coords[2 * i], coords[2 * i + 1]);
        if (d < minDist && d > 0) {
            i1 = i;
            minDist = d;
        }
    }

    var minRadius = Infinity;

    // find the third point which forms the smallest circumcircle with the first two
    for (i = 0; i < points.length; i++) {
        if (i === i0 || i === i1) continue;

        var r = circumradius(
            coords[2 * i0], coords[2 * i0 + 1],
            coords[2 * i1], coords[2 * i1 + 1],
            coords[2 * i], coords[2 * i + 1]);

        if (r < minRadius) {
            i2 = i;
            minRadius = r;
        }
    }

    if (minRadius === Infinity) {
        throw new Error('No Delaunay triangulation exists for this input.');
    }

    // swap the order of the seed points for counter-clockwise orientation
    if (area(coords[2 * i0], coords[2 * i0 + 1],
        coords[2 * i1], coords[2 * i1 + 1],
        coords[2 * i2], coords[2 * i2 + 1]) < 0) {

        var tmp = i1;
        i1 = i2;
        i2 = tmp;
    }

    var i0x = coords[2 * i0];
    var i0y = coords[2 * i0 + 1];
    var i1x = coords[2 * i1];
    var i1y = coords[2 * i1 + 1];
    var i2x = coords[2 * i2];
    var i2y = coords[2 * i2 + 1];

    var center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
    this._cx = center.x;
    this._cy = center.y;

    // sort the points by distance from the seed triangle circumcenter
    quicksort(ids, coords, 0, ids.length - 1, center.x, center.y);

    // initialize a hash table for storing edges of the advancing convex hull
    this._hashSize = Math.ceil(Math.sqrt(points.length));
    this._hash = [];
    for (i = 0; i < this._hashSize; i++) this._hash[i] = null;

    // initialize a circular doubly-linked list that will hold an advancing convex hull
    var e = this.hull = insertNode(coords, i0);
    this._hashEdge(e);
    e.t = 0;
    e = insertNode(coords, i1, e);
    this._hashEdge(e);
    e.t = 1;
    e = insertNode(coords, i2, e);
    this._hashEdge(e);
    e.t = 2;

    var maxTriangles = 2 * points.length - 5;
    var triangles = this.triangles = new Uint32Array(maxTriangles * 3);
    var halfedges = this.halfedges = new Int32Array(maxTriangles * 3);

    this.trianglesLen = 0;

    this._addTriangle(i0, i1, i2, -1, -1, -1);

    var xp, yp;
    for (var k = 0; k < ids.length; k++) {
        i = ids[k];
        x = coords[2 * i];
        y = coords[2 * i + 1];

        // skip duplicate points
        if (x === xp && y === yp) continue;
        xp = x;
        yp = y;

        // skip seed triangle points
        if ((x === i0x && y === i0y) ||
            (x === i1x && y === i1y) ||
            (x === i2x && y === i2y)) continue;

        // find a visible edge on the convex hull using edge hash
        var startKey = this._hashKey(x, y);
        var key = startKey;
        var start;
        do {
            start = this._hash[key];
            key = (key + 1) % this._hashSize;
        } while ((!start || start.removed) && key !== startKey);

        e = start;
        while (area(x, y, e.x, e.y, e.next.x, e.next.y) >= 0) {
            e = e.next;
            if (e === start) {
                throw new Error('Something is wrong with the input points.');
            }
        }

        var walkBack = e === start;

        // add the first triangle from the point
        var t = this._addTriangle(e.i, i, e.next.i, -1, -1, e.t);

        e.t = t; // keep track of boundary triangles on the hull
        e = insertNode(coords, i, e);

        // recursively flip triangles from the point until they satisfy the Delaunay condition
        e.t = this._legalize(t + 2);
        if (e.prev.prev.t === halfedges[t + 1]) {
            e.prev.prev.t = t + 2;
        }

        // walk forward through the hull, adding more triangles and flipping recursively
        var q = e.next;
        while (area(x, y, q.x, q.y, q.next.x, q.next.y) < 0) {
            t = this._addTriangle(q.i, i, q.next.i, q.prev.t, -1, q.t);
            q.prev.t = this._legalize(t + 2);
            this.hull = removeNode(q);
            q = q.next;
        }

        if (walkBack) {
            // walk backward from the other side, adding more triangles and flipping
            q = e.prev;
            while (area(x, y, q.prev.x, q.prev.y, q.x, q.y) < 0) {
                t = this._addTriangle(q.prev.i, i, q.i, -1, q.t, q.prev.t);
                this._legalize(t + 2);
                q.prev.t = t;
                this.hull = removeNode(q);
                q = q.prev;
            }
        }

        // save the two new edges in the hash table
        this._hashEdge(e);
        this._hashEdge(e.prev);
    }

    // trim typed triangle mesh arrays
    this.triangles = triangles.subarray(0, this.trianglesLen);
    this.halfedges = halfedges.subarray(0, this.trianglesLen);
}

Delaunator.prototype = {

    _hashEdge: function (e) {
        this._hash[this._hashKey(e.x, e.y)] = e;
    },

    _hashKey: function (x, y) {
        var dx = x - this._cx;
        var dy = y - this._cy;
        // use pseudo-angle: a measure that monotonically increases
        // with real angle, but doesn't require expensive trigonometry
        var p = 1 - dx / (Math.abs(dx) + Math.abs(dy));
        return Math.floor((2 + (dy < 0 ? -p : p)) / 4 * this._hashSize);
    },

    _legalize: function (a) {
        var triangles = this.triangles;
        var coords = this.coords;
        var halfedges = this.halfedges;

        var b = halfedges[a];

        var a0 = a - a % 3;
        var b0 = b - b % 3;

        var al = a0 + (a + 1) % 3;
        var ar = a0 + (a + 2) % 3;
        var br = b0 + (b + 1) % 3;
        var bl = b0 + (b + 2) % 3;

        var p0 = triangles[ar];
        var pr = triangles[a];
        var pl = triangles[al];
        var p1 = triangles[bl];

        var illegal = inCircle(
            coords[2 * p0], coords[2 * p0 + 1],
            coords[2 * pr], coords[2 * pr + 1],
            coords[2 * pl], coords[2 * pl + 1],
            coords[2 * p1], coords[2 * p1 + 1]);

        if (illegal) {
            triangles[a] = p1;
            triangles[b] = p0;

            this._link(a, halfedges[bl]);
            this._link(b, halfedges[ar]);
            this._link(ar, bl);

            this._legalize(a);
            return this._legalize(br);
        }

        return ar;
    },

    _link: function (a, b) {
        this.halfedges[a] = b;
        if (b !== -1) this.halfedges[b] = a;
    },

    // add a new triangle given vertex indices and adjacent half-edge ids
    _addTriangle: function (i0, i1, i2, a, b, c) {
        var t = this.trianglesLen;

        this.triangles[t] = i0;
        this.triangles[t + 1] = i1;
        this.triangles[t + 2] = i2;

        this._link(t, a);
        this._link(t + 1, b);
        this._link(t + 2, c);

        this.trianglesLen += 3;

        return t;
    }
};

function dist(ax, ay, bx, by) {
    var dx = ax - bx;
    var dy = ay - by;
    return dx * dx + dy * dy;
}

function area(px, py, qx, qy, rx, ry) {
    return (qy - py) * (rx - qx) - (qx - px) * (ry - qy);
}

function inCircle(ax, ay, bx, by, cx, cy, px, py) {
    ax -= px;
    ay -= py;
    bx -= px;
    by -= py;
    cx -= px;
    cy -= py;

    var ap = ax * ax + ay * ay;
    var bp = bx * bx + by * by;
    var cp = cx * cx + cy * cy;

    return ax * (by * cp - bp * cy) -
           ay * (bx * cp - bp * cx) +
           ap * (bx * cy - by * cx) < 0;
}

function circumradius(ax, ay, bx, by, cx, cy) {
    bx -= ax;
    by -= ay;
    cx -= ax;
    cy -= ay;

    var bl = bx * bx + by * by;
    var cl = cx * cx + cy * cy;

    if (bl === 0 || cl === 0) return Infinity;

    var d = bx * cy - by * cx;
    if (d === 0) return Infinity;

    var x = (cy * bl - by * cl) * 0.5 / d;
    var y = (bx * cl - cx * bl) * 0.5 / d;

    return x * x + y * y;
}

function circumcenter(ax, ay, bx, by, cx, cy) {
    bx -= ax;
    by -= ay;
    cx -= ax;
    cy -= ay;

    var bl = bx * bx + by * by;
    var cl = cx * cx + cy * cy;

    var d = bx * cy - by * cx;

    var x = (cy * bl - by * cl) * 0.5 / d;
    var y = (bx * cl - cx * bl) * 0.5 / d;

    return {
        x: ax + x,
        y: ay + y
    };
}

// create a new node in a doubly linked list
function insertNode(coords, i, prev) {
    var node = {
        i: i,
        x: coords[2 * i],
        y: coords[2 * i + 1],
        t: 0,
        prev: null,
        next: null,
        removed: false
    };

    if (!prev) {
        node.prev = node;
        node.next = node;

    } else {
        node.next = prev.next;
        node.prev = prev;
        prev.next.prev = node;
        prev.next = node;
    }
    return node;
}

function removeNode(node) {
    node.prev.next = node.next;
    node.next.prev = node.prev;
    node.removed = true;
    return node.prev;
}

function quicksort(ids, coords, left, right, cx, cy) {
    var i, j, temp;

    if (right - left <= 20) {
        for (i = left + 1; i <= right; i++) {
            temp = ids[i];
            j = i - 1;
            while (j >= left && compare(coords, ids[j], temp, cx, cy) > 0) ids[j + 1] = ids[j--];
            ids[j + 1] = temp;
        }
    } else {
        var median = (left + right) >> 1;
        i = left + 1;
        j = right;
        swap(ids, median, i);
        if (compare(coords, ids[left], ids[right], cx, cy) > 0) swap(ids, left, right);
        if (compare(coords, ids[i], ids[right], cx, cy) > 0) swap(ids, i, right);
        if (compare(coords, ids[left], ids[i], cx, cy) > 0) swap(ids, left, i);

        temp = ids[i];
        while (true) {
            do i++; while (compare(coords, ids[i], temp, cx, cy) < 0);
            do j--; while (compare(coords, ids[j], temp, cx, cy) > 0);
            if (j < i) break;
            swap(ids, i, j);
        }
        ids[left + 1] = ids[j];
        ids[j] = temp;

        if (right - i + 1 >= j - left) {
            quicksort(ids, coords, i, right, cx, cy);
            quicksort(ids, coords, left, j - 1, cx, cy);
        } else {
            quicksort(ids, coords, left, j - 1, cx, cy);
            quicksort(ids, coords, i, right, cx, cy);
        }
    }
}

function compare(coords, i, j, cx, cy) {
    var d1 = dist(coords[2 * i], coords[2 * i + 1], cx, cy);
    var d2 = dist(coords[2 * j], coords[2 * j + 1], cx, cy);
    return (d1 - d2) || (coords[2 * i] - coords[2 * j]) || (coords[2 * i + 1] - coords[2 * j + 1]);
}

function swap(arr, i, j) {
    var tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

function defaultGetX(p) {
    return p[0];
}
function defaultGetY(p) {
    return p[1];
}

},{}],4:[function(require,module,exports){
"use strict"

function iota(n) {
  var result = new Array(n)
  for(var i=0; i<n; ++i) {
    result[i] = i
  }
  return result
}

module.exports = iota
},{}],5:[function(require,module,exports){
/*!
 * Determine if an object is a Buffer
 *
 * @author   Feross Aboukhadijeh <feross@feross.org> <http://feross.org>
 * @license  MIT
 */

// The _isBuffer check is for Safari 5-7 support, because it's missing
// Object.prototype.constructor. Remove this eventually
module.exports = function (obj) {
  return obj != null && (isBuffer(obj) || isSlowBuffer(obj) || !!obj._isBuffer)
}

function isBuffer (obj) {
  return !!obj.constructor && typeof obj.constructor.isBuffer === 'function' && obj.constructor.isBuffer(obj)
}

// For Node v0.10 support. Remove this eventually.
function isSlowBuffer (obj) {
  return typeof obj.readFloatLE === 'function' && typeof obj.slice === 'function' && isBuffer(obj.slice(0, 0))
}

},{}],6:[function(require,module,exports){
module.exports = moore

function moore(range, dims) {
  dims = dims || 2
  range = range || 1
  return recurse([], [], 0)

  function recurse(array, temp, d) {
    if (d === dims-1) {
      for (var i = -range; i <= range; i += 1) {
        if (i || temp.some(function(n) {
          return n
        })) array.push(temp.concat(i))
      }
    } else {
      for (var i = -range; i <= range; i += 1) {
        recurse(array, temp.concat(i), d+1)
      }
    }
    return array
  }
}

},{}],7:[function(require,module,exports){
var iota = require("iota-array")
var isBuffer = require("is-buffer")

var hasTypedArrays  = ((typeof Float64Array) !== "undefined")

function compare1st(a, b) {
  return a[0] - b[0]
}

function order() {
  var stride = this.stride
  var terms = new Array(stride.length)
  var i
  for(i=0; i<terms.length; ++i) {
    terms[i] = [Math.abs(stride[i]), i]
  }
  terms.sort(compare1st)
  var result = new Array(terms.length)
  for(i=0; i<result.length; ++i) {
    result[i] = terms[i][1]
  }
  return result
}

function compileConstructor(dtype, dimension) {
  var className = ["View", dimension, "d", dtype].join("")
  if(dimension < 0) {
    className = "View_Nil" + dtype
  }
  var useGetters = (dtype === "generic")

  if(dimension === -1) {
    //Special case for trivial arrays
    var code =
      "function "+className+"(a){this.data=a;};\
var proto="+className+".prototype;\
proto.dtype='"+dtype+"';\
proto.index=function(){return -1};\
proto.size=0;\
proto.dimension=-1;\
proto.shape=proto.stride=proto.order=[];\
proto.lo=proto.hi=proto.transpose=proto.step=\
function(){return new "+className+"(this.data);};\
proto.get=proto.set=function(){};\
proto.pick=function(){return null};\
return function construct_"+className+"(a){return new "+className+"(a);}"
    var procedure = new Function(code)
    return procedure()
  } else if(dimension === 0) {
    //Special case for 0d arrays
    var code =
      "function "+className+"(a,d) {\
this.data = a;\
this.offset = d\
};\
var proto="+className+".prototype;\
proto.dtype='"+dtype+"';\
proto.index=function(){return this.offset};\
proto.dimension=0;\
proto.size=1;\
proto.shape=\
proto.stride=\
proto.order=[];\
proto.lo=\
proto.hi=\
proto.transpose=\
proto.step=function "+className+"_copy() {\
return new "+className+"(this.data,this.offset)\
};\
proto.pick=function "+className+"_pick(){\
return TrivialArray(this.data);\
};\
proto.valueOf=proto.get=function "+className+"_get(){\
return "+(useGetters ? "this.data.get(this.offset)" : "this.data[this.offset]")+
"};\
proto.set=function "+className+"_set(v){\
return "+(useGetters ? "this.data.set(this.offset,v)" : "this.data[this.offset]=v")+"\
};\
return function construct_"+className+"(a,b,c,d){return new "+className+"(a,d)}"
    var procedure = new Function("TrivialArray", code)
    return procedure(CACHED_CONSTRUCTORS[dtype][0])
  }

  var code = ["'use strict'"]

  //Create constructor for view
  var indices = iota(dimension)
  var args = indices.map(function(i) { return "i"+i })
  var index_str = "this.offset+" + indices.map(function(i) {
        return "this.stride[" + i + "]*i" + i
      }).join("+")
  var shapeArg = indices.map(function(i) {
      return "b"+i
    }).join(",")
  var strideArg = indices.map(function(i) {
      return "c"+i
    }).join(",")
  code.push(
    "function "+className+"(a," + shapeArg + "," + strideArg + ",d){this.data=a",
      "this.shape=[" + shapeArg + "]",
      "this.stride=[" + strideArg + "]",
      "this.offset=d|0}",
    "var proto="+className+".prototype",
    "proto.dtype='"+dtype+"'",
    "proto.dimension="+dimension)

  //view.size:
  code.push("Object.defineProperty(proto,'size',{get:function "+className+"_size(){\
return "+indices.map(function(i) { return "this.shape["+i+"]" }).join("*"),
"}})")

  //view.order:
  if(dimension === 1) {
    code.push("proto.order=[0]")
  } else {
    code.push("Object.defineProperty(proto,'order',{get:")
    if(dimension < 4) {
      code.push("function "+className+"_order(){")
      if(dimension === 2) {
        code.push("return (Math.abs(this.stride[0])>Math.abs(this.stride[1]))?[1,0]:[0,1]}})")
      } else if(dimension === 3) {
        code.push(
"var s0=Math.abs(this.stride[0]),s1=Math.abs(this.stride[1]),s2=Math.abs(this.stride[2]);\
if(s0>s1){\
if(s1>s2){\
return [2,1,0];\
}else if(s0>s2){\
return [1,2,0];\
}else{\
return [1,0,2];\
}\
}else if(s0>s2){\
return [2,0,1];\
}else if(s2>s1){\
return [0,1,2];\
}else{\
return [0,2,1];\
}}})")
      }
    } else {
      code.push("ORDER})")
    }
  }

  //view.set(i0, ..., v):
  code.push(
"proto.set=function "+className+"_set("+args.join(",")+",v){")
  if(useGetters) {
    code.push("return this.data.set("+index_str+",v)}")
  } else {
    code.push("return this.data["+index_str+"]=v}")
  }

  //view.get(i0, ...):
  code.push("proto.get=function "+className+"_get("+args.join(",")+"){")
  if(useGetters) {
    code.push("return this.data.get("+index_str+")}")
  } else {
    code.push("return this.data["+index_str+"]}")
  }

  //view.index:
  code.push(
    "proto.index=function "+className+"_index(", args.join(), "){return "+index_str+"}")

  //view.hi():
  code.push("proto.hi=function "+className+"_hi("+args.join(",")+"){return new "+className+"(this.data,"+
    indices.map(function(i) {
      return ["(typeof i",i,"!=='number'||i",i,"<0)?this.shape[", i, "]:i", i,"|0"].join("")
    }).join(",")+","+
    indices.map(function(i) {
      return "this.stride["+i + "]"
    }).join(",")+",this.offset)}")

  //view.lo():
  var a_vars = indices.map(function(i) { return "a"+i+"=this.shape["+i+"]" })
  var c_vars = indices.map(function(i) { return "c"+i+"=this.stride["+i+"]" })
  code.push("proto.lo=function "+className+"_lo("+args.join(",")+"){var b=this.offset,d=0,"+a_vars.join(",")+","+c_vars.join(","))
  for(var i=0; i<dimension; ++i) {
    code.push(
"if(typeof i"+i+"==='number'&&i"+i+">=0){\
d=i"+i+"|0;\
b+=c"+i+"*d;\
a"+i+"-=d}")
  }
  code.push("return new "+className+"(this.data,"+
    indices.map(function(i) {
      return "a"+i
    }).join(",")+","+
    indices.map(function(i) {
      return "c"+i
    }).join(",")+",b)}")

  //view.step():
  code.push("proto.step=function "+className+"_step("+args.join(",")+"){var "+
    indices.map(function(i) {
      return "a"+i+"=this.shape["+i+"]"
    }).join(",")+","+
    indices.map(function(i) {
      return "b"+i+"=this.stride["+i+"]"
    }).join(",")+",c=this.offset,d=0,ceil=Math.ceil")
  for(var i=0; i<dimension; ++i) {
    code.push(
"if(typeof i"+i+"==='number'){\
d=i"+i+"|0;\
if(d<0){\
c+=b"+i+"*(a"+i+"-1);\
a"+i+"=ceil(-a"+i+"/d)\
}else{\
a"+i+"=ceil(a"+i+"/d)\
}\
b"+i+"*=d\
}")
  }
  code.push("return new "+className+"(this.data,"+
    indices.map(function(i) {
      return "a" + i
    }).join(",")+","+
    indices.map(function(i) {
      return "b" + i
    }).join(",")+",c)}")

  //view.transpose():
  var tShape = new Array(dimension)
  var tStride = new Array(dimension)
  for(var i=0; i<dimension; ++i) {
    tShape[i] = "a[i"+i+"]"
    tStride[i] = "b[i"+i+"]"
  }
  code.push("proto.transpose=function "+className+"_transpose("+args+"){"+
    args.map(function(n,idx) { return n + "=(" + n + "===undefined?" + idx + ":" + n + "|0)"}).join(";"),
    "var a=this.shape,b=this.stride;return new "+className+"(this.data,"+tShape.join(",")+","+tStride.join(",")+",this.offset)}")

  //view.pick():
  code.push("proto.pick=function "+className+"_pick("+args+"){var a=[],b=[],c=this.offset")
  for(var i=0; i<dimension; ++i) {
    code.push("if(typeof i"+i+"==='number'&&i"+i+">=0){c=(c+this.stride["+i+"]*i"+i+")|0}else{a.push(this.shape["+i+"]);b.push(this.stride["+i+"])}")
  }
  code.push("var ctor=CTOR_LIST[a.length+1];return ctor(this.data,a,b,c)}")

  //Add return statement
  code.push("return function construct_"+className+"(data,shape,stride,offset){return new "+className+"(data,"+
    indices.map(function(i) {
      return "shape["+i+"]"
    }).join(",")+","+
    indices.map(function(i) {
      return "stride["+i+"]"
    }).join(",")+",offset)}")

  //Compile procedure
  var procedure = new Function("CTOR_LIST", "ORDER", code.join("\n"))
  return procedure(CACHED_CONSTRUCTORS[dtype], order)
}

function arrayDType(data) {
  if(isBuffer(data)) {
    return "buffer"
  }
  if(hasTypedArrays) {
    switch(Object.prototype.toString.call(data)) {
      case "[object Float64Array]":
        return "float64"
      case "[object Float32Array]":
        return "float32"
      case "[object Int8Array]":
        return "int8"
      case "[object Int16Array]":
        return "int16"
      case "[object Int32Array]":
        return "int32"
      case "[object Uint8Array]":
        return "uint8"
      case "[object Uint16Array]":
        return "uint16"
      case "[object Uint32Array]":
        return "uint32"
      case "[object Uint8ClampedArray]":
        return "uint8_clamped"
    }
  }
  if(Array.isArray(data)) {
    return "array"
  }
  return "generic"
}

var CACHED_CONSTRUCTORS = {
  "float32":[],
  "float64":[],
  "int8":[],
  "int16":[],
  "int32":[],
  "uint8":[],
  "uint16":[],
  "uint32":[],
  "array":[],
  "uint8_clamped":[],
  "buffer":[],
  "generic":[]
}

;(function() {
  for(var id in CACHED_CONSTRUCTORS) {
    CACHED_CONSTRUCTORS[id].push(compileConstructor(id, -1))
  }
});

function wrappedNDArrayCtor(data, shape, stride, offset) {
  if(data === undefined) {
    var ctor = CACHED_CONSTRUCTORS.array[0]
    return ctor([])
  } else if(typeof data === "number") {
    data = [data]
  }
  if(shape === undefined) {
    shape = [ data.length ]
  }
  var d = shape.length
  if(stride === undefined) {
    stride = new Array(d)
    for(var i=d-1, sz=1; i>=0; --i) {
      stride[i] = sz
      sz *= shape[i]
    }
  }
  if(offset === undefined) {
    offset = 0
    for(var i=0; i<d; ++i) {
      if(stride[i] < 0) {
        offset -= (shape[i]-1)*stride[i]
      }
    }
  }
  var dtype = arrayDType(data)
  var ctor_list = CACHED_CONSTRUCTORS[dtype]
  while(ctor_list.length <= d+1) {
    ctor_list.push(compileConstructor(dtype, ctor_list.length-1))
  }
  var ctor = ctor_list[d+1]
  return ctor(data, shape, stride, offset)
}

module.exports = wrappedNDArrayCtor

},{"iota-array":4,"is-buffer":5}],8:[function(require,module,exports){
"use strict";

module.exports = require('./src/poisson-disk-sampling');

},{"./src/poisson-disk-sampling":10}],9:[function(require,module,exports){
"use strict";

module.exports = function euclideanDistanceN (point1, point2) {
    var result = 0,
        i = 0;

    for (; i < point1.length; i++) {
        result += Math.pow(point1[i] - point2[i], 2);
    }

    return Math.sqrt(result);
};

},{}],10:[function(require,module,exports){
"use strict";

var zeros = require('zeros'),
    moore = require('moore'),
    euclideanDistanceN = require('./euclidean-distance'),
    sphereRandom = require('./sphere-random');

/**
 * Get the neighbourhood ordered by distance, including the origin point
 * @param {int} dimensionNumber Number of dimensions
 * @returns {Array} Neighbourhood
 */
var getNeighbourhood = function getNeighbourhood (dimensionNumber) {
    var neighbourhood = moore(2, dimensionNumber),
        origin = [],
        dimension;

    for (dimension = 0; dimension < dimensionNumber; dimension++) {
        origin.push(0);
    }

    neighbourhood.push(origin);

    // sort by ascending distance to optimize proximity checks
    // see point 5.1 in Parallel Poisson Disk Sampling by Li-Yi Wei, 2008
    // http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.460.3061&rank=1
    neighbourhood.sort(function (n1, n2) {
        var squareDist1 = 0,
            squareDist2 = 0;

        for (var dimension = 0; dimension < dimensionNumber; dimension++) {
            squareDist1 += Math.pow(n1[dimension], 2);
            squareDist2 += Math.pow(n2[dimension], 2);
        }

        if (squareDist1 < squareDist2) {
            return -1;
        } else if(squareDist1 > squareDist2) {
            return 1;
        } else {
            return 0;
        }
    });

    return neighbourhood;
};


/**
 * PoissonDiskSampling constructor
 * @param {Array} shape Shape of the space
 * @param {float} minDistance Minimum distance between each points
 * @param {float} [maxDistance] Maximum distance between each points, defaults to minDistance * 2
 * @param {int} [maxTries] Number of times the algorithm has to try to place a point in the neighbourhood of another points, defaults to 30
 * @param {function|null} [rng] RNG function, defaults to Math.random
 * @constructor
 */
var PoissonDiskSampling = function PoissonDiskSampling (shape, minDistance, maxDistance, maxTries, rng) {
    maxDistance = maxDistance || minDistance * 2;

    this.shape = shape;
    this.dimension = this.shape.length;
    this.minDistance = minDistance;
    this.deltaDistance = maxDistance - minDistance;
    this.cellSize = minDistance / Math.sqrt(this.dimension);
    this.maxTries = maxTries || 30;
    this.rng = rng || Math.random;

    this.neighbourhood = getNeighbourhood(this.dimension);

    this.currentPoint = null;
    this.processList = [];
    this.samplePoints = [];

    // cache grid

    this.gridShape = [];

    for (var i = 0; i < this.dimension; i++) {
        this.gridShape.push(Math.ceil(shape[i] / this.cellSize));
    }

    this.grid = zeros(this.gridShape, 'uint32'); //will store references to samplePoints
};

PoissonDiskSampling.prototype.shape = null;
PoissonDiskSampling.prototype.dimension = null;
PoissonDiskSampling.prototype.minDistance = null;
PoissonDiskSampling.prototype.deltaDistance = null;
PoissonDiskSampling.prototype.cellSize = null;
PoissonDiskSampling.prototype.maxTries = null;
PoissonDiskSampling.prototype.rng = null;
PoissonDiskSampling.prototype.neighbourhood = null;

PoissonDiskSampling.prototype.currentPoint = null;
PoissonDiskSampling.prototype.processList = null;
PoissonDiskSampling.prototype.samplePoints = null;
PoissonDiskSampling.prototype.gridShape = null;
PoissonDiskSampling.prototype.grid = null;

/**
 * Add a totally random point in the grid
 * @returns {Array} The point added to the grid
 */
PoissonDiskSampling.prototype.addRandomPoint = function () {
    var point = new Array(this.dimension);

    for (var i = 0; i < this.dimension; i++) {
        point[i] = this.rng() * this.shape[i];
    }

    return this.directAddPoint(point);
};

/**
 * Add a given point to the grid
 * @param {Array} point Point
 * @returns {Array|null} The point added to the grid, null if the point is out of the bound or not of the correct dimension
 */
PoissonDiskSampling.prototype.addPoint = function (point) {
    var dimension,
        valid = true;

    if (point.length === this.dimension) {
        for (dimension = 0; dimension < this.dimension && valid; dimension++) {
            valid = (point[dimension] >= 0 && point[dimension] <= this.shape[dimension]);
        }
    } else {
        valid = false;
    }

    return valid ? this.directAddPoint(point) : null;
};

/**
 * Add a given point to the grid, without any check
 * @param {Array} point Point
 * @returns {Array} The point added to the grid
 * @protected
 */
PoissonDiskSampling.prototype.directAddPoint = function (point) {
    var internalArrayIndex = 0,
        stride = this.grid.stride,
        dimension;

    this.processList.push(point);
    this.samplePoints.push(point);

    for (dimension = 0; dimension < this.dimension; dimension++) {
        internalArrayIndex += ((point[dimension] / this.cellSize) | 0) * stride[dimension];
    }

    this.grid.data[internalArrayIndex] = this.samplePoints.length; // store the point reference

    return point;
};

/**
 * Check whether a given point is in the neighbourhood of existing points
 * @param {Array} point Point
 * @returns {boolean} Whether the point is in the neighbourhood of another point
 * @protected
 */
PoissonDiskSampling.prototype.inNeighbourhood = function (point) {
    var dimensionNumber = this.dimension,
        stride = this.grid.stride,
        neighbourIndex,
        internalArrayIndex,
        dimension,
        currentDimensionValue,
        existingPoint;

    for (neighbourIndex = 0; neighbourIndex < this.neighbourhood.length; neighbourIndex++) {
        internalArrayIndex = 0;

        for (dimension = 0; dimension < dimensionNumber; dimension++) {
            currentDimensionValue = ((point[dimension] / this.cellSize) | 0) + this.neighbourhood[neighbourIndex][dimension];

            if (currentDimensionValue >= 0 && currentDimensionValue < this.gridShape[dimension]) {
                internalArrayIndex += currentDimensionValue * stride[dimension];
            }
        }

        if (this.grid.data[internalArrayIndex] !== 0) {
            existingPoint = this.samplePoints[this.grid.data[internalArrayIndex] - 1];

            if (euclideanDistanceN(point, existingPoint) < this.minDistance) {
                return true;
            }
        }
    }

    return false;
};

/**
 * Try to generate a new point in the grid, returns null if it wasn't possible
 * @returns {Array|null} The added point or null
 */
PoissonDiskSampling.prototype.next = function () {
    var tries,
        angle,
        distance,
        currentPoint,
        newPoint,
        inShape,
        i;

    while (this.processList.length > 0) {
        if (this.currentPoint === null) {
            this.currentPoint = this.processList.shift();
        }

        currentPoint = this.currentPoint;

        for (tries = 0; tries < this.maxTries; tries++) {
            inShape = true;
            distance = this.minDistance + this.deltaDistance * this.rng();

            if (this.dimension === 2) {
                angle = this.rng() * Math.PI * 2;
                newPoint = [
                    Math.cos(angle),
                    Math.sin(angle)
                ];
            } else {
                newPoint = sphereRandom(this.dimension, this.rng);
            }

            for (i = 0; inShape && i < this.dimension; i++) {
                newPoint[i] = currentPoint[i] + newPoint[i] * distance;
                inShape = (newPoint[i] >= 0 && newPoint[i] <= this.shape[i] - 1)
            }

            if (inShape && !this.inNeighbourhood(newPoint)) {
                return this.directAddPoint(newPoint);
            }
        }

        if (tries >= this.maxTries) {
            this.currentPoint = null;
        }
    }

    return null;
};

/**
 * Automatically fill the grid, adding a random point to start the process if needed.
 * Will block the thread, probably best to use it in a web worker or child process.
 * @returns {Array[]} Sample points
 */
PoissonDiskSampling.prototype.fill = function () {
    if (this.samplePoints.length === 0) {
        this.addRandomPoint();
    }

    while(this.next()) {}

    return this.samplePoints;
};

/**
 * Get all the points in the grid.
 * @returns {Array[]} Sample points
 */
PoissonDiskSampling.prototype.getAllPoints = function () {
    return this.samplePoints;
};

/**
 * Reinitialize the grid as well as the internal state
 */
PoissonDiskSampling.prototype.reset = function () {
    var gridData = this.grid.data,
        i = 0;

    // reset the cache grid
    for (i = 0; i < gridData.length; i++) {
        gridData[i] = 0;
    }

    // new array for the samplePoints as it is passed by reference to the outside
    this.samplePoints = [];

    // reset the internal state
    this.currentPoint = null;
    this.processList.length = 0;
};

module.exports = PoissonDiskSampling;

},{"./euclidean-distance":9,"./sphere-random":11,"moore":6,"zeros":12}],11:[function(require,module,exports){
"use strict";

// sphere-random module by Mikola Lysenko under the MIT License
// waiting for https://github.com/scijs/sphere-random/pull/1 to be merged

module.exports = sampleSphere;

var defaultRng = Math.random;

/**
 * @param {int} d Dimensions
 * @param {Function} [rng]
 * @returns {Array}
 */
function sampleSphere(d, rng) {
    var v = new Array(d),
        d2 = Math.floor(d/2) << 1,
        r2 = 0.0,
        rr,
        r,
        theta,
        h,
        i;

    rng = rng || defaultRng;

    for (i = 0; i < d2; i += 2) {
        rr = -2.0 * Math.log(rng());
        r =  Math.sqrt(rr);
        theta = 2.0 * Math.PI * rng();

        r2+= rr;
        v[i] = r * Math.cos(theta);
        v[i+1] = r * Math.sin(theta);
    }

    if (d % 2) {
        var x = Math.sqrt(-2.0 * Math.log(rng())) * Math.cos(2.0 * Math.PI * rng());
        v[d - 1] = x;
        r2+= Math.pow(x, 2);
    }

    h = 1.0 / Math.sqrt(r2);

    for (i=0; i<d; ++i) {
        v[i] *= h;
    }

    return v;
}

},{}],12:[function(require,module,exports){
"use strict"

var ndarray = require("ndarray")

function dtypeToType(dtype) {
  switch(dtype) {
    case 'uint8':
      return Uint8Array;
    case 'uint16':
      return Uint16Array;
    case 'uint32':
      return Uint32Array;
    case 'int8':
      return Int8Array;
    case 'int16':
      return Int16Array;
    case 'int32':
      return Int32Array;
    case 'float':
    case 'float32':
      return Float32Array;
    case 'double':
    case 'float64':
      return Float64Array;
    case 'uint8_clamped':
      return Uint8ClampedArray;
    case 'generic':
    case 'buffer':
    case 'data':
    case 'dataview':
      return ArrayBuffer;
    case 'array':
      return Array;
  }
}

module.exports = function zeros(shape, dtype) {
  dtype = dtype || 'float64';
  var sz = 1;
  for(var i=0; i<shape.length; ++i) {
    sz *= shape[i];
  }
  return ndarray(new (dtypeToType(dtype))(sz), shape);
}

},{"ndarray":7}],13:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */

'use strict';

const util = require('./util');

function biome(ocean, water, coast, temperature, moisture) {
    if (ocean) {
        return 'OCEAN';
    } else if (water) {
        if (temperature > 0.9) return 'MARSH';
        if (temperature < 0.2) return 'ICE';
        return 'LAKE';
    } else if (coast) {
        return 'BEACH';
    } else if (temperature < 0.2) {
        if (moisture > 0.50) return 'SNOW';
        else if (moisture > 0.33) return 'TUNDRA';
        else if (moisture > 0.16) return 'BARE';
        else return 'SCORCHED';
    } else if (temperature < 0.4) {
        if (moisture > 0.66) return 'TAIGA';
        else if (moisture > 0.33) return 'SHRUBLAND';
        else return 'TEMPERATE_DESERT';
    } else if (temperature < 0.7) {
        if (moisture > 0.83) return 'TEMPERATE_RAIN_FOREST';
        else if (moisture > 0.50) return 'TEMPERATE_DECIDUOUS_FOREST';
        else if (moisture > 0.16) return 'GRASSLAND';
        else return 'TEMPERATE_DESERT';
    } else {
        if (moisture > 0.66) return 'TROPICAL_RAIN_FOREST';
        else if (moisture > 0.33) return 'TROPICAL_SEASONAL_FOREST';
        else if (moisture > 0.16) return 'GRASSLAND';
        else return 'SUBTROPICAL_DESERT';
    }
}


/**
 * A coast region is land that has an ocean neighbor
 */
exports.assign_r_coast = function(r_coast, mesh, r_ocean) {
    r_coast.length = mesh.numRegions;
    r_coast.fill(false);
    
    let out_r = [];
    for (let r1 = 0; r1 < mesh.numRegions; r1++) {
        mesh.r_circulate_r(out_r, r1);
        if (!r_ocean[r1]) {
            for (let r2 of out_r) {
                if (r_ocean[r2]) {
                    r_coast[r1] = true;
                    break;
                }
            }
        }
    }
    return r_coast;
};


/**
 * Temperature assignment
 *
 * Temperature is based on elevation and latitude.
 * The normal range is 0.0=cold, 1.0=hot, but it is not 
 * limited to that range, especially when using temperature bias.
 *
 * The northernmost parts of the map get bias_north added to them;
 * the southernmost get bias_south added; in between it's a blend.
 */
exports.assign_r_temperature = function(
    r_temperature,
    mesh,
    r_ocean, r_water,
    r_elevation, r_moisture,
    bias_north, bias_south
) {
    r_temperature.length = mesh.numRegions;
    for (let r = 0; r < mesh.numRegions; r++) {
        let latitude = mesh.r_y(r) / 1000; /* 0.0 - 1.0 */
        let d_temperature = util.mix(bias_north, bias_south, latitude);
        r_temperature[r] = 1.0 - r_elevation[r] + d_temperature;
    }
    return r_temperature;
};


/**
 * Biomes assignment -- see the biome() function above
 */
exports.assign_r_biome = function(
    r_biome,
    mesh,
    r_ocean, r_water, r_coast, r_temperature, r_moisture
) {
    r_biome.length = mesh.numRegions;
    for (let r = 0; r < mesh.numRegions; r++) {
        r_biome[r] = biome(r_ocean[r], r_water[r], r_coast[r],
                           r_temperature[r], r_moisture[r]);
    }
    return r_biome;
};

},{"./util":19}],14:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */

'use strict';

/**
 * Coast corners are connected to coast sides, which have
 * ocean on one side and land on the other
 */
function find_coasts_t(mesh, r_ocean) {
    let coasts_t = [];
    for (let s = 0; s < mesh.numSides; s++) {
        let r0 = mesh.s_begin_r(s);
        let r1 = mesh.s_end_r(s);
        let t = mesh.s_inner_t(s);
        if (r_ocean[r0] && !r_ocean[r1]) {
            // It might seem that we also need to check !r_ocean[r0] && r_ocean[r1]
            // and it might seem that we have to add both t and its opposite but
            // each t vertex shows up in *four* directed sides, so we only have to test
            // one fourth of those conditions to get the vertex in the list once.
            coasts_t.push(t);
        }
    }
    return coasts_t;
}


/**
 * Elevation is based on breadth first search from the seed points,
 * which are the coastal graph nodes. Since breadth first search also
 * calculates the 'parent' pointers, return those for use as the downslope
 * graph. To handle lakes, which should have all corners at the same elevation,
 * there are two deviations from breadth first search:
 * 1. Instead of pushing to the end of the queue, push to the beginning.
 * 2. Like uniform cost search, check if the new distance is better than
 *    previously calculated distances. It is possible that one lake corner
 *    was reached with distance 2 and another with distance 3, and we need
 *    to revisit that node and make sure it's set to 2.
 */
exports.assign_t_elevation = function(
    t_elevation, t_coastdistance, t_downslope_s,
    mesh,
    r_ocean, r_water, randInt
) {
    t_coastdistance.length = mesh.numTriangles;
    t_downslope_s.length = mesh.numTriangles;
    t_elevation.length = mesh.numTriangles;
    t_coastdistance.fill(null);
    t_downslope_s.fill(-1);
    
    const t_ocean = (t) => r_ocean[mesh.s_begin_r(3*t)];
    const r_lake = (r) => r_water[r] && !r_ocean[r];
    const s_lake = (s) => r_lake(mesh.s_begin_r(s)) || r_lake(mesh.s_end_r(s));

    let out_s = [];
    let queue_t = find_coasts_t(mesh, r_ocean);
    queue_t.forEach((t) => { t_coastdistance[t] = 0; });
    let minDistance = 1, maxDistance = 1;
    
    while (queue_t.length > 0) {
        let current_t = queue_t.shift();
        mesh.t_circulate_s(out_s, current_t);
        let iOffset = randInt(out_s.length);
        for (let i = 0; i < out_s.length; i++) {
            let s = out_s[(i + iOffset) % out_s.length];
            let lake = s_lake(s);
            let neighbor_t = mesh.s_outer_t(s);
            let newDistance = (lake? 0 : 1) + t_coastdistance[current_t];
            if (t_coastdistance[neighbor_t] === null || newDistance < t_coastdistance[neighbor_t]) {
                t_downslope_s[neighbor_t] = mesh.s_opposite_s(s);
                t_coastdistance[neighbor_t] = newDistance;
                if (t_ocean(neighbor_t) && newDistance > minDistance) { minDistance = newDistance; }
                if (!t_ocean(neighbor_t) && newDistance > maxDistance) { maxDistance = newDistance; }
                if (lake) {
                    queue_t.unshift(neighbor_t);
                } else {
                    queue_t.push(neighbor_t);
                }
            }
        }
    }

    t_coastdistance.forEach((d, t) => {
        t_elevation[t] = t_ocean(t) ? (-d / minDistance) : (d / maxDistance);
    });
};


/** 
 * Set r elevation to the average of the t elevations. There's a
 * corner case though: it is possible for an ocean region (r) to be
 * surrounded by coastline corners (t), and coastlines are set to 0
 * elevation. This means the region elevation would be 0. To avoid
 * this, I subtract a small amount for ocean regions. */
exports.assign_r_elevation = function(r_elevation, mesh, t_elevation, r_ocean) {
    const max_ocean_elevation = -0.01;
    r_elevation.length = mesh.numRegions;
    let out_t = [];
    for (let r = 0; r < mesh.numRegions; r++) {
        mesh.r_circulate_t(out_t, r);
        let elevation = 0.0;
        for (let t of out_t) {
            elevation += t_elevation[t];
        }
        r_elevation[r] = elevation/out_t.length;
        if (r_ocean[r] && r_elevation[r] > max_ocean_elevation) {
            r_elevation[r] = max_ocean_elevation;
        }
    }
    return r_elevation;
};


/**
 * Redistribute elevation values so that lower elevations are more common
 * than higher elevations. Specifically, we want elevation Z to have frequency
 * (1-Z), for all the non-ocean regions.
 */
// TODO: this messes up lakes, as they will no longer all be at the same elevation
exports.redistribute_t_elevation = function(t_elevation, mesh) {
    // NOTE: This is the same algorithm I used in 2010, because I'm
    // trying to recreate that map generator to some extent. I don't
    // think it's a great approach for other games but it worked well
    // enough for that one.
    
    // SCALE_FACTOR increases the mountain area. At 1.0 the maximum
    // elevation barely shows up on the map, so we set it to 1.1.
    const SCALE_FACTOR = 1.1;

    let nonocean_t = [];
    for (let t = 0; t < mesh.numSolidTriangles; t++) {
        if (t_elevation[t] > 0.0) {
            nonocean_t.push(t);
        }
    }
    
    nonocean_t.sort((t1, t2) => t_elevation[t1] - t_elevation[t2]);

    for (let i = 0; i < nonocean_t.length; i++) {
        // Let y(x) be the total area that we want at elevation <= x.
        // We want the higher elevations to occur less than lower
        // ones, and set the area to be y(x) = 1 - (1-x)^2.
        let y = i / (nonocean_t.length-1);
        // Now we have to solve for x, given the known y.
        //  *  y = 1 - (1-x)^2
        //  *  y = 1 - (1 - 2x + x^2)
        //  *  y = 2x - x^2
        //  *  x^2 - 2x + y = 0
        // From this we can use the quadratic equation to get:
        let x = Math.sqrt(SCALE_FACTOR) - Math.sqrt(SCALE_FACTOR*(1-y));
        if (x > 1.0) x = 1.0;
        t_elevation[nonocean_t[i]] = x;
    }
};

},{}],15:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *      http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

const util =         require('./util');
const Water =        require('./water');
const Elevation =    require('./elevation');
const Rivers =       require('./rivers');
const Moisture =     require('./moisture');
const Biomes =       require('./biomes');
const NoisyEdges =   require('./noisy-edges');

/**
 * Map generator
 *
 * Map coordinates are 0 ≤ x ≤ 1000, 0 ≤ y ≤ 1000.
 *
 * mesh: DualMesh
 * noisyEdgeOptions: {length, amplitude, seed}
 * makeRandInt: function(seed) -> function(N) -> an int from 0 to N-1
 */
class Map {
    constructor(mesh, noisyEdgeOptions, makeRandInt) {
        this.mesh = mesh;
        this.makeRandInt = makeRandInt;
        this.s_lines = NoisyEdges.assign_s_segments(
            [],
            this.mesh,
            noisyEdgeOptions,
            this.makeRandInt(noisyEdgeOptions.seed)
        );

        this.r_water = [];
        this.r_ocean = [];
        this.t_coastdistance = [];
        this.t_elevation = [];
        this.t_downslope_s = [];
        this.r_elevation = [];
        this.s_flow = [];
        this.r_waterdistance = [];
        this.r_moisture = [];
        this.r_coast = [];
        this.r_temperature = [];
        this.r_biome = [];
    }

 
    calculate(options) {
        options = Object.assign({
            noise: null, // required: function(nx, ny) -> number from -1 to +1
            shape: {round: 0.5, inflate: 0.4, amplitudes: [1/2, 1/4, 1/8, 1/16]},
            numRivers: 30,
            drainageSeed: 0,
            riverSeed: 0,
            noisyEdge: {length: 10, amplitude: 0.2, seed: 0},
            biomeBias: {north_temperature: 0, south_temperature: 0, moisture: 0},
        }, options);

        Water.assign_r_water(this.r_water, this.mesh, options.noise, options.shape);
        Water.assign_r_ocean(this.r_ocean, this.mesh, this.r_water);
        
        Elevation.assign_t_elevation(
            this.t_elevation, this.t_coastdistance, this.t_downslope_s,
            this.mesh,
            this.r_ocean, this.r_water, this.makeRandInt(options.drainageSeed)
        );
        Elevation.redistribute_t_elevation(this.t_elevation, this.mesh);
        Elevation.assign_r_elevation(this.r_elevation, this.mesh, this.t_elevation, this.r_ocean);

        this.spring_t = Rivers.find_spring_t(this.mesh, this.r_water, this.t_elevation, this.t_downslope_s);
        util.randomShuffle(this.spring_t, this.makeRandInt(options.riverSeed));
        
        this.river_t = this.spring_t.slice(0, options.numRivers);
        Rivers.assign_s_flow(this.s_flow, this.mesh, this.t_downslope_s, this.river_t, this.t_elevation);
        
        Moisture.assign_r_moisture(
            this.r_moisture, this.r_waterdistance,
            this.mesh,
            this.r_water, Moisture.find_moisture_seeds_r(this.mesh, this.s_flow, this.r_ocean, this.r_water)
        );
        Moisture.redistribute_r_moisture(this.r_moisture, this.mesh, this.r_water,
                                         options.biomeBias.moisture, 1 + options.biomeBias.moisture);

        Biomes.assign_r_coast(this.r_coast, this.mesh, this.r_ocean);
        Biomes.assign_r_temperature(
            this.r_temperature,
            this.mesh,
            this.r_ocean, this.r_water, this.r_elevation, this.r_moisture,
            options.biomeBias.north_temperature, options.biomeBias.south_temperature
        );
        Biomes.assign_r_biome(
            this.r_biome,
            this.mesh,
            this.r_ocean, this.r_water, this.r_coast, this.r_temperature, this.r_moisture
        );
    }
}

module.exports = Map;

},{"./biomes":13,"./elevation":14,"./moisture":16,"./noisy-edges":17,"./rivers":18,"./util":19,"./water":20}],16:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */

'use strict';

/**
 * Find regions adjacent to rivers; out_r should be a Set
 */
exports.find_riverbanks_r = function(out_r, mesh, s_flow) {
    for (let s = 0; s < mesh.numSolidSides; s++) {
        if (s_flow[s] > 0) {
            out_r.add(mesh.s_begin_r(s));
            out_r.add(mesh.s_end_r(s));
        }
    }
};


/**
 * Find lakeshores -- regions adjacent to lakes; out_r should be a Set
 */
exports.find_lakeshores_r = function(out_r, mesh, r_ocean, r_water) {
    for (let s = 0; s < mesh.numSolidSides; s++) {
        let r0 = mesh.s_begin_r(s),
            r1 = mesh.s_end_r(s);
        if (r_water[r0] && !r_ocean[r0]) {
            out_r.add(r0);
            out_r.add(r1);
        }
    }
};


/**
 * Find regions that have maximum moisture; returns a Set
 */
exports.find_moisture_seeds_r = function(mesh, s_flow, r_ocean, r_water) {
    let seeds_r = new Set();
    exports.find_riverbanks_r(seeds_r, mesh, s_flow);
    exports.find_lakeshores_r(seeds_r, mesh, r_ocean, r_water);
    return seeds_r;
};


/**
 * Assign moisture level. Oceans and lakes have moisture 1.0. Land
 * regions have moisture based on the distance to the nearest fresh
 * water. Lakeshores and riverbanks are distance 0. Moisture will be
 * 1.0 at distance 0 and go down to 0.0 at the maximum distance.
 */
exports.assign_r_moisture = function(
    r_moisture, r_waterdistance,
    mesh,
    r_water, seed_r /* Set */
) {
    r_waterdistance.length = mesh.numRegions;
    r_moisture.length = mesh.numRegions;
    r_waterdistance.fill(null);
    
    let out_r = [];
    let queue_r = Array.from(seed_r);
    let maxDistance = 1;
    queue_r.forEach((r) => { r_waterdistance[r] = 0; });
    while (queue_r.length > 0) {
        let current_r = queue_r.shift();
        mesh.r_circulate_r(out_r, current_r);
        for (let neighbor_r of out_r) {
            if (!r_water[neighbor_r] && r_waterdistance[neighbor_r] === null) {
                let newDistance = 1 + r_waterdistance[current_r];
                r_waterdistance[neighbor_r] = newDistance;
                if (newDistance > maxDistance) { maxDistance = newDistance; }
                queue_r.push(neighbor_r);
            }
        }
    }

    r_waterdistance.forEach((d, r) => {
        r_moisture[r] = r_water[r]? 1.0 : 1.0 - Math.pow(d / maxDistance, 0.5);
    });
};


/**
 * Redistribute moisture values evenly so that all moistures
 * from min_moisture to max_moisture are equally represented.
 */
exports.redistribute_r_moisture = function(r_moisture, mesh, r_water, min_moisture, max_moisture) {
    let land_r = [];
    for (let r = 0; r < mesh.numSolidRegions; r++) {
        if (!r_water[r]) {
            land_r.push(r);
        }
    }

    land_r.sort((r1, r2) => r_moisture[r1] - r_moisture[r2]);
    
    for (let i = 0; i < land_r.length; i++) {
        r_moisture[land_r[i]] = min_moisture + (max_moisture-min_moisture) * i / (land_r.length - 1);
    }
};

},{}],17:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */

'use strict';

const {mixp} = require('./util');

/**
 * Noisy edges is a variant of midpoint subdivision that keeps the lines
 * constrained to a quadrilateral. See the explanation here:
 * http://www.redblobgames.com/maps/mapgen2/noisy-edges.html
 */

/**
 * Return the noisy line from a to b, within quadrilateral a-p-b-q,
 * as an array of points, not including a. The recursive subdivision
 * has up to 2^levels segments. Segments below a given length are
 * not subdivided further.
 */
const divisor = 0x10000000;
exports.recursiveSubdivision = (length, amplitude, randInt) =>
    function recur(a, b, p, q) {
        let dx = a[0] - b[0], dy = a[1] - b[1];
        if (dx*dx + dy*dy < length*length) { return [b]; }
        
        let ap = mixp([], a, p, 0.5),
            bp = mixp([], b, p, 0.5),
            aq = mixp([], a, q, 0.5),
            bq = mixp([], b, q, 0.5);

        let division = 0.5 * (1 - amplitude) + randInt(divisor)/divisor * amplitude;
        let center = mixp([], p, q, division);
        
        let results1 = recur(a, center, ap, aq),
            results2 = recur(center, b, bp, bq);

        return results1.concat(results2);
    };


// TODO: this allocates lots of tiny arrays; find a data format that
// doesn't have so many allocations

exports.assign_s_segments = function(s_lines, mesh, {amplitude, length}, randInt) {
    s_lines.length = mesh.numSides;
    for (let s = 0; s < mesh.numSides; s++) {
        let t0 = mesh.s_inner_t(s),
            t1 = mesh.s_outer_t(s),
            r0 = mesh.s_begin_r(s),
            r1 = mesh.s_end_r(s);
        if (r0 < r1) {
            if (mesh.s_ghost(s)) {
                s_lines[s] = [mesh.t_pos([], t1)];
            } else {
                s_lines[s] = exports.recursiveSubdivision(length, amplitude, randInt)(
                    mesh.t_pos([], t0),
                    mesh.t_pos([], t1),
                    mesh.r_pos([], r0),
                    mesh.r_pos([], r1)
                );
            }
            // construct line going the other way; since the line is a
            // half-open interval with [p1, p2, p3, ..., pn] but not
            // p0, we want to reverse all but the last element, and
            // then append p0
            let opposite = s_lines[s].slice(0, -1);
            opposite.reverse();
            opposite.push(mesh.t_pos([], t0));
            s_lines[mesh.s_opposite_s(s)] = opposite;
        }
    }
    return s_lines;
};

},{"./util":19}],18:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */

'use strict';

const MIN_SPRING_ELEVATION = 0.3;
const MAX_SPRING_ELEVATION = 0.9;

/**
 * Find candidates for river sources
 *
 * Unlike the assign_* functions this does not write into an existing array
 */
exports.find_spring_t = function(mesh, r_water, t_elevation, t_downslope_s) {
    const t_water = (t) =>
          (  r_water[mesh.s_begin_r(3*t)]
          || r_water[mesh.s_begin_r(3*t+1)]
          || r_water[mesh.s_begin_r(3*t+2)] );

    let spring_t = new Set();
    // Add everything above some elevation, but not lakes
    for (let t = 0; t < mesh.numSolidTriangles; t++) {
        if (t_elevation[t] >= MIN_SPRING_ELEVATION &&
            t_elevation[t] <= MAX_SPRING_ELEVATION &&
            !t_water(t)) {
            spring_t.add(t);
        }
    }
    return Array.from(spring_t);
};


exports.assign_s_flow = function(s_flow, mesh, t_downslope_s, river_t) {
    // Each river in river_t contributes 1 flow down to the coastline
    s_flow.length = mesh.numSides;
    s_flow.fill(0);
    for (let t of river_t) {
        for (;;) {
            let s = t_downslope_s[t];
            if (s === -1) { break; }
            s_flow[s]++;
            let next_t = mesh.s_outer_t(s);
            if (next_t === t) { break; }
            t = next_t;
        }
    }
    return s_flow;
};

},{}],19:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */

/**
 * Return value, unless it's undefined, then return orElse 
 */
exports.fallback = function(value, orElse) {
    return (value !== undefined)? value : orElse;
};

/**
 * Add several noise values together
 */
exports.fbm_noise = function(noise, amplitudes, nx, ny) {
    let sum = 0, sumOfAmplitudes = 0;
    for (let octave = 0; octave < amplitudes.length; octave++) {
        let frequency = 1 << octave;
        sum += amplitudes[octave] * noise.noise2D(nx * frequency, ny * frequency, octave);
        sumOfAmplitudes += amplitudes[octave];
    }
    return sum / sumOfAmplitudes;
};

/**
 * Like GLSL. Return t clamped to the range [lo,hi] inclusive 
 */
exports.clamp = function(t, lo, hi) {
    if (t < lo) { return lo; }
    if (t > hi) { return hi; }
    return t;
};

/**
 * Like GLSL. Return a mix of a and b; all a when is 0 and all b when
 * t is 1; extrapolates when t outside the range [0,1] 
 */
exports.mix = function(a, b, t) {
    return a * (1.0-t) + b * t;
};

/**
 * Componentwise mix for arrays of equal length; output goes in 'out'
 */
exports.mixp = function(out, p, q, t) {
    out.length = p.length;
    for (let i = 0; i < p.length; i++) {
        out[i] = exports.mix(p[i], q[i], t);
    }
    return out;
};

/**
 * Like GLSL. 
 */
exports.smoothstep = function(a, b, t) {
    // https://en.wikipedia.org/wiki/Smoothstep
    if (t <= a) { return 0; }
    if (t >= b) { return 1; }
    t = (t - a) / (b - a);
    return (3 - 2*t) * t * t;
};

/**
 * Circumcenter of a triangle with vertices a,b,c
 */
exports.circumcenter = function(a, b, c) {
    // https://en.wikipedia.org/wiki/Circumscribed_circle#Circumcenter_coordinates
    let ad = a[0]*a[0] + a[1]*a[1],
        bd = b[0]*b[0] + b[1]*b[1],
        cd = c[0]*c[0] + c[1]*c[1];
    let D = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    let Ux = 1/D * (ad * (b[1] - c[1]) + bd * (c[1] - a[1]) + cd * (a[1] - b[1]));
    let Uy = 1/D * (ad * (c[0] - b[0]) + bd * (a[0] - c[0]) + cd * (b[0] - a[0]));
    return [Ux, Uy];
};

/**
 * Intersection of line p1--p2 and line p3--p4,
 * between 0.0 and 1.0 if it's in the line segment
 */
exports.lineIntersection = function(x1, y1, x2, y2, x3, y3, x4, y4) {
    // from http://paulbourke.net/geometry/pointlineplane/
    let ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));
    let ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));
    return {ua, ub};
};

/**
 * in-place shuffle of an array - Fisher-Yates
 * https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle#The_modern_algorithm
 */
exports.randomShuffle = function(array, randInt) {
    for (let i = array.length-1; i > 0; i--) {
        let j = randInt(i+1);
        let swap = array[i];
        array[i] = array[j];
        array[j] = swap;
    }
    return array;
};

},{}],20:[function(require,module,exports){
// From http://www.redblobgames.com/maps/mapgen2/
// Copyright 2017 Red Blob Games <redblobgames@gmail.com>
// License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>

'use strict';
const util = require('./util');

// NOTE: r_water, r_ocean, other fields are boolean valued so it
// could be more efficient to pack them as bit fields in Uint8Array

/* a region is water if the noise value is low */
exports.assign_r_water = function(r_water, mesh, noise, params) {
    r_water.length = mesh.numRegions;
    for (let r = 0; r < mesh.numRegions; r++) {
        if (mesh.r_ghost(r) || mesh.r_boundary(r)) {
            r_water[r] = true;
        } else {
            let nx = (mesh.r_x(r) - 500) / 500;
            let ny = (mesh.r_y(r) - 500) / 500;
            let distance = Math.max(Math.abs(nx), Math.abs(ny));
            let n = util.fbm_noise(noise, params.amplitudes, nx, ny);
            n = util.mix(n, 0.5, params.round);
            r_water[r] = n - (1.0 - params.inflate) * distance*distance < 0;
        }
    }
    return r_water;
};


/* a region is ocean if it is a water region connected to the ghost region,
   which is outside the boundary of the map; this could be any seed set but
   for islands, the ghost region is a good seed */
exports.assign_r_ocean = function(r_ocean, mesh, r_water) {
    r_ocean.length = mesh.numRegions;
    r_ocean.fill(false);
    let stack = [mesh.ghost_r()];
    let r_out = [];
    while (stack.length > 0) {
        let r1 = stack.pop();
        mesh.r_circulate_r(r_out, r1);
        for (let r2 of r_out) {
            if (r_water[r2] && !r_ocean[r2]) {
                r_ocean[r2] = true;
                stack.push(r2);
            }
        }
    }
    return r_ocean;
};

},{"./util":19}],21:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */
'use strict';

const hashInt = require('hash-int');

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

},{"hash-int":22}],22:[function(require,module,exports){
"use strict"

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
}

module.exports = hashInt

},{}],23:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */
'use strict';

exports.discreteColors = {
    OCEAN: "#44447a",
    COAST: "#33335a",
    LAKESHORE: "#225588",
    LAKE: "#336699",
    RIVER: "#225588",
    MARSH: "#2f6666",
    ICE: "#99ffff",
    BEACH: "#a09077",
    SNOW: "#ffffff",
    TUNDRA: "#bbbbaa",
    BARE: "#888888",
    SCORCHED: "#555555",
    TAIGA: "#99aa77",
    SHRUBLAND: "#889977",
    TEMPERATE_DESERT: "#c9d29b",
    TEMPERATE_RAIN_FOREST: "#448855",
    TEMPERATE_DECIDUOUS_FOREST: "#679459",
    GRASSLAND: "#88aa55",
    SUBTROPICAL_DESERT: "#d2b98b",
    TROPICAL_RAIN_FOREST: "#337755",
    TROPICAL_SEASONAL_FOREST: "#559944",
};

function smoothColoring(e, t, m) {
    // adapted from <https://www.redblobgames.com/maps/terrain-from-noise/>
    if (e < 0.0) {
        return `rgb(${(48 + 48*e) | 0}, ${(64 + 64*e) | 0}, ${(127 + 128*e) | 0})`;
    }

    // Green or brown at low elevation, and make it more white-ish
    // as you get colder
    let white = (1-t) * (1-t);
    m = 1.0 - ((1-m)*(1-m));
    var red = 210 - 100*m, grn = 185 - 45*m, blu = 139 - 45*m;
    return `rgb(${(255 * white + red * (1-white)) | 0}, 
                ${(255 * white + grn * (1-white)) | 0}, 
                ${(255 * white + blu * (1-white)) | 0})`;
}


class Coloring {
    constructor() {
    }

    draw_coast_s(map, s) {
        return map.r_ocean[map.mesh.s_begin_r(s)] !== map.r_ocean[map.mesh.s_end_r(s)];
    }

    draw_lakeside_s(map, s) {
        let r0 = map.mesh.s_begin_r(s),
            r1 = map.mesh.s_end_r(s);
        return (map.r_water[r0] !== map.r_water[r1]
                && !map.r_ocean[r0]
                && map.r_biome[r0] !== 'ICE'
                && map.r_biome[r1] !== 'ICE');
    }
    
    draw_river_s(map, s) {
        let r0 = map.mesh.s_begin_r(s),
            r1 = map.mesh.s_end_r(s);
        return ((map.s_flow[s] > 0 || map.s_flow[map.mesh.s_opposite_s(s)] > 0)
                && !map.r_water[r0] && !map.r_water[r1]);
    }

    biome(map, r) {
        return "red";
    }

    side(map, s) {
        let r0 = map.mesh.s_begin_r(s),
            r1 = map.mesh.s_end_r(s);
        if (this.draw_coast_s(map, s)) {
            // Coastlines are thick
            return {
                noisy: true,
                lineWidth: 3,
                strokeStyle: exports.discreteColors.COAST,
            };
        } else if (this.draw_lakeside_s(map, s)) {
            // Lake boundary
            return {
                noisy: true,
                lineWidth: 1.5,
                strokeStyle: exports.discreteColors.LAKESHORE,
            };
        } else if (this.draw_river_s(map, s)) {
            // River
            return {
                noisy: true,
                lineWidth: 2.0 * Math.sqrt(map.s_flow[s]),
                strokeStyle: exports.discreteColors.RIVER,
            };
        } else if (map.r_biome[r0] === map.r_biome[r1]) {
            return {
                noisy: false,
                lineWidth: 1.0,
                strokeStyle: this.biome(map, r0),
            };
        } else {
            return {
                noisy: true,
                lineWidth: 1.0,
                strokeStyle: this.biome(map, r0),
            };
        }
    }
}

class Discrete extends Coloring {
    biome(map, r) {
        return exports.discreteColors[map.r_biome[r]];
    }
}

class Smooth extends Coloring {
    biome(map, r) {
        if (map.r_water[r] && !map.r_ocean[r]) {
            return exports.discreteColors[map.r_biome[r]];
        } else {
            return smoothColoring(
                map.r_elevation[r],
                Math.min(1, Math.max(0, map.r_temperature[r])),
                Math.min(1, Math.max(0, map.r_moisture[r]))
            );
        }
    }
}

exports.Discrete = Discrete;
exports.Smooth = Smooth;

},{}],24:[function(require,module,exports){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 * License: Apache v2.0 <http://www.apache.org/licenses/LICENSE-2.0.html>
 */
'use strict';

const util = require('@redblobgames/mapgen2/util');
const Colormap = require('./colormap');

const noiseSize = 100;
let noiseCanvas = null;
function makeNoise(randInt) {
    if (noiseCanvas === null) {
        noiseCanvas = document.createElement('canvas');
        noiseCanvas.width = noiseSize;
        noiseCanvas.height = noiseSize;

        let ctx = noiseCanvas.getContext('2d');
        const imageData = ctx.getImageData(0, 0, noiseSize, noiseSize);
        const pixels = imageData.data;

        for (let y = 0, p = 0; y < noiseSize; y++) {
            for (let x = 0; x < noiseSize; x++) {
                let value = 128 + randInt(16) - 8;
                pixels[p++] = value;
                pixels[p++] = value;
                pixels[p++] = value;
                pixels[p++] = 255;
            }
        }
        ctx.putImageData(imageData, 0, 0);
    }
}


exports.noisyFill = function(ctx, width, height, randInt) {
    makeNoise(randInt);
    ctx.save();
    ctx.globalCompositeOperation = 'soft-light';
    ctx.drawImage(noiseCanvas, 0, 0, width, height);
    ctx.globalCompositeOperation = 'hard-light';
    for (let y = 0; y < height; y += noiseSize) {
        for (let x = 0; x < width; x += noiseSize) {
            ctx.drawImage(noiseCanvas, x, y, noiseSize, noiseSize);
        }
    }
    ctx.restore();
};


const lightSize = 250;
const lightScaleZ = 15;
const lightVector = [-1, -1, 0];
let lightCanvas = null;

// quick & dirty light based on normal vector
function calculateLight(ax, ay, az,
                        bx, by, bz,
                        cx, cy, cz) {
    az *= lightScaleZ;
    bz *= lightScaleZ;
    cz *= lightScaleZ;
    let ux = bx - ax, uy = by - ay, uz = bz - az,
        vx = cx - ax, vy = cy - ay, vz = cz - az;
    // cross product (ugh I should have a lib for this)
    let nx = uy*vz - uz*vy,
        ny = uz*vx - ux*vz,
        nz = ux*vy - uy*vx;
    let length = -Math.sqrt(nx*nx + ny*ny + nz*nz);
    nx /= length;
    ny /= length;
    nz /= length;
    let dotProduct = nx * lightVector[0] + ny * lightVector[1] + nz * lightVector[2];
    let light = 0.5 + 10 * dotProduct;
    return util.clamp(light, 0, 1);
}

function makeLight(map) {
    if (lightCanvas === null) {
        lightCanvas = document.createElement('canvas');
        lightCanvas.width = lightSize;
        lightCanvas.height = lightSize;
    }
    let ctx = lightCanvas.getContext('2d');
    ctx.save();
    ctx.scale(lightSize/1000, lightSize/1000);
    ctx.fillStyle = "hsl(0,0%,50%)";
    ctx.fillRect(0, 0, 1000, 1000);
    let mesh = map.mesh;

    // Draw lighting on land; skip in the ocean
    let r_out = [];
    for (let t = 0; t < mesh.numSolidTriangles; t++) {
        mesh.t_circulate_r(r_out, t);
        if (r_out.some((r) => map.r_water[r])) { continue; }
        let ax = mesh.r_x(r_out[0]),
            ay = mesh.r_y(r_out[0]),
            az = map.r_elevation[r_out[0]],
            bx = mesh.r_x(r_out[1]),
            by = mesh.r_y(r_out[1]),
            bz = map.r_elevation[r_out[1]],
            cx = mesh.r_x(r_out[2]),
            cy = mesh.r_y(r_out[2]),
            cz = map.r_elevation[r_out[2]];
        let light = calculateLight(ax, ay, az*az, bx, by, bz*bz, cx, cy, cz*cz);
        light = util.mix(light, map.t_elevation[t], 0.5);
        ctx.strokeStyle = ctx.fillStyle = `hsl(0,0%,${(light*100) | 0}%)`;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(ax, ay);
        ctx.lineTo(bx, by);
        ctx.lineTo(cx, cy);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
    }

    ctx.restore();
}
        
exports.lighting = function(ctx, width, height, map) {
    makeLight(map);
    ctx.globalCompositeOperation = 'soft-light';
    ctx.drawImage(lightCanvas, 0, 0, width, height);
}


const islandShapeSize = 200;
let islandShapeCanvas = null;
function makeIsland(noise, params) {
    if (!islandShapeCanvas) {
        islandShapeCanvas = document.createElement('canvas');
        islandShapeCanvas.width = islandShapeSize;
        islandShapeCanvas.height = islandShapeSize;
    }
    
    let ctx = islandShapeCanvas.getContext('2d');
    const imageData = ctx.getImageData(0, 0, islandShapeSize, islandShapeSize);
    const pixels = imageData.data;

    for (let y = 0, p = 0; y < islandShapeSize; y++) {
        let ny = 2 * y/islandShapeSize - 1;
        for (let x = 0; x < islandShapeSize; x++) {
            let nx = 2 * x/islandShapeSize - 1;
            let distance = Math.max(Math.abs(nx), Math.abs(ny));
            let n = util.fbm_noise(noise, params.amplitudes, nx, ny);
            n = util.mix(n, 0.5, params.round);
            if (n - (1.0 - params.inflate) * distance*distance < 0) {
                // water color uses OCEAN discrete color
                pixels[p++] = 0x44;
                pixels[p++] = 0x44;
                pixels[p++] = 0x7a;
            } else {
                // land color uses BEACH discrete color
                pixels[p++] = 0xa0;
                pixels[p++] = 0x90;
                pixels[p++] = 0x77;
            }
            pixels[p++] = 255;
        }
    }
    
    ctx.putImageData(imageData, 0, 0);
}

exports.approximateIslandShape = function(ctx, width, height, noise, params) {
    makeIsland(noise, params);
    ctx.drawImage(islandShapeCanvas, 0, 0, width, height);
};


exports.background = function(ctx, colormap) {
    ctx.fillStyle = Colormap.discreteColors.OCEAN;
    ctx.fillRect(0, 0, 1000, 1000);
};


exports.noisyRegions = function(ctx, map, colormap, noisyEdge) {
    let {mesh} = map;
    let out_s = [];

    for (let r = 0; r < mesh.numSolidRegions; r++) {
        mesh.r_circulate_s(out_s, r);
        let last_t = mesh.s_inner_t(out_s[0]);
        ctx.fillStyle = ctx.strokeStyle = colormap.biome(map, r);
        ctx.beginPath();
        ctx.moveTo(mesh.t_x(last_t), mesh.t_y(last_t));
        for (let s of out_s) {
            if (!noisyEdge || !colormap.side(map, s).noisy) {
                let first_t = mesh.s_outer_t(s);
                ctx.lineTo(mesh.t_x(first_t), mesh.t_y(first_t));
            } else {
                for (let p of map.s_lines[s]) {
                    ctx.lineTo(p[0], p[1]);
                }
            }
        }
        ctx.fill();
    }
};

/*
 * Helper function: how big is the region?
 *
 * Returns the minimum distance from the region center to a corner
 */
function region_radius(mesh, r) {
    let rx = mesh.r_x(r), ry = mesh.r_y(r);
    let min_distance_squared = Infinity;
    let out_t = [];
    mesh.r_circulate_t(out_t, r);
    for (let t of out_t) {
        let tx = mesh.t_x(t), ty = mesh.t_y(t);
        let dx = rx - tx, dy = ry - ty;
        let distance_squared = dx*dx + dy*dy;
        if (distance_squared < min_distance_squared) {
            min_distance_squared = distance_squared;
        }
    }
    return Math.sqrt(min_distance_squared);
}

/*
 * Draw a biome icon in each of the regions
 */
exports.regionIcons = function(ctx, map, mapIconsConfig, randInt) {
    let {mesh} = map;
    for (let r = 0; r < mesh.numSolidRegions; r++) {
        if (mesh.r_boundary(r)) { continue; }
        let biome = map.r_biome[r];
        let radius = region_radius(mesh, r);
        let row = {
            OCEAN: 0, LAKE: 0,
            SHRUBLAND: 2,
            TEMPERATE_DESERT: 3, SUBTROPICAL_DESERT: 3,
            TROPICAL_RAIN_FOREST: 4, TROPICAL_SEASONAL_FOREST: 4,
            TEMPERATE_DECIDUOUS_FOREST: 5, TEMPERATE_RAIN_FOREST: 5,
            GRASSLAND: 6,
            MARSH: 7,
            TAIGA: 9,
        }[biome];
        // NOTE: mountains reflect elevation, but the biome
        // calculation reflects temperature, so if you set the biome
        // bias to be 'cold', you'll get more snow, but you shouldn't
        // get more mountains, so the mountains are calculated
        // separately from biomes
        if (row === 5 && mesh.r_y(r) < 300) { row = 9; }
        if (map.r_elevation[r] > 0.8) { row = 1; }
        if (row === undefined) { continue; }
        let col = 1 + randInt(5);
        ctx.drawImage(mapIconsConfig.image,
                      mapIconsConfig.left + col*100, mapIconsConfig.top + row*100,
                      100, 100,
                      mesh.r_x(r) - radius, mesh.r_y(r) - radius,
                      2*radius, 2*radius);
    }
};


/*
 * Drawing the region polygons leaves little gaps in HTML5 Canvas
 * so I need to draw edges to fill those gaps. Sometimes those edges
 * are simple straight lines but sometimes they're thick noisy lines
 * like coastlines and rivers.
 *
 * This step is rather slow so it's split up into phases.
 *
 * If 'filter' is defined, filter(side, style) should return true if
 * the edge is to be drawn. This is used by the rivers and coastline
 * drawing functions.
 */
exports.noisyEdges = function(ctx, map, colormap, noisyEdge, phase /* 0-15 */, filter=null) {
    let {mesh} = map;
    let begin = (mesh.numSolidSides/16 * phase) | 0;
    let end = (mesh.numSolidSides/16 * (phase+1)) | 0;
    for (let s = begin; s < end; s++) {
        let style = colormap.side(map, s);
        if (filter && !filter(s, style)) { continue; }
        ctx.strokeStyle = style.strokeStyle;
        ctx.lineWidth = style.lineWidth;
        let last_t = mesh.s_inner_t(s);
        ctx.beginPath();
        ctx.moveTo(mesh.t_x(last_t), mesh.t_y(last_t));
        if (!noisyEdge || !style.noisy) {
            let first_t = mesh.s_outer_t(s);
            ctx.lineTo(mesh.t_x(first_t), mesh.t_y(first_t));
        } else {
            for (let p of map.s_lines[s]) {
                ctx.lineTo(p[0], p[1]);
            }
        }
        ctx.stroke();
    }
};


exports.vertices = function(ctx, map) {
    let {mesh} = map;
    ctx.fillStyle = "black";
    for (let r = 0; r < mesh.numSolidRegions; r++) {
        ctx.beginPath();
        ctx.arc(mesh.r_x(r), mesh.r_y(r), 2, 0, 2*Math.PI);
        ctx.fill();
    }
};


exports.rivers = function(ctx, map, colormap, noisyEdge, fast) {
    if (!fast) {
        ctx.lineCap = 'round';
        ctx.lineJoin = 'round';
    }
    for (let phase = 0; phase < 16; phase++) {
        exports.noisyEdges(ctx, map, colormap, noisyEdge, phase,
                           (s, style) => colormap.draw_river_s(map, s));
    }
};


exports.coastlines = function(ctx, map, colormap, noisyEdge) {
    ctx.lineCap = 'round';
    ctx.lineJoin = 'round';
    for (let phase = 0; phase < 16; phase++) {
        exports.noisyEdges(ctx, map, colormap, noisyEdge, phase,
                           (s, style) => colormap.draw_coast_s(map, s));
    }
};

},{"./colormap":23,"@redblobgames/mapgen2/util":19}],25:[function(require,module,exports){
(function (global){
/*
 * From http://www.redblobgames.com/maps/mapgen2/
 * Copyright 2017 Red Blob Games <redblobgames@gmail.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *      http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
'use strict';

const SimplexNoise = require('simplex-noise');
const DualMesh =     require('@redblobgames/dual-mesh');
const createMesh =   require('@redblobgames/dual-mesh/create');
const Map =          require('@redblobgames/mapgen2');
const Draw =         require('./draw');
const Colormap =     require('./colormap');
const urlUtils =     require('url-search-utils');
const {makeRandInt, makeRandFloat} = require('@redblobgames/prng');


let defaultUiState = {
    seed: 187,
    variant: 0,
    size: 'medium',
    'noisy-fills': true,
    'noisy-edges': true,
    icons: true,
    biomes: false,
    lighting: false,
    'north-temperature': 0,
    'south-temperature': 0,
    rainfall: 0,
    canvasSize: 0,
    persistence: 0, /* 0 means "normal" persistence value of 1/2 */
};

let uiState = {};
Object.assign(uiState, defaultUiState);


let _mapCache = [];
function getMap(size) {
    const spacing = {
        tiny: 38,
        small: 26,
        medium: 18,
        large: 12.8,
        huge: 9,
    };
    if (!_mapCache[size]) {
        // NOTE: the seeds here are constant so that I can reuse the same
        // mesh and noisy edges for all maps, but you could get more variety
        // by creating a new Map object each time
        _mapCache[size] = new Map(
            new DualMesh(createMesh({spacing: spacing[size], random: makeRandFloat(12345)})),
            {amplitude: 0.2, length: 4, seed: 12345},
            makeRandInt
        );
        console.log(`Map size "${size}" has ${_mapCache[size].mesh.numRegions} regions`);
    }
    return _mapCache[size];
}


/**
 * Manage drawing with requestAnimationFrame
 *
 * 1. Each frame call one function from the queue.
 * 2. If the queue empties, stop calling requestAnimationFrame.
 */
const processingPerFrameInMs = 1000/60;
let requestAnimationFrameId = null;
let requestAnimationFrameQueue = [];
function requestAnimationFrameHandler() {
    requestAnimationFrameId = null;
    let timeStart = performance.now();
    while (requestAnimationFrameQueue.length > 0
           && performance.now() - timeStart < processingPerFrameInMs) {
        let f = requestAnimationFrameQueue.shift();
        f();
    }
    if (requestAnimationFrameQueue.length > 0) {
        requestAnimationFrameId = requestAnimationFrame(requestAnimationFrameHandler);
    }
}


/* map icons */
const mapIconsConfig = {left: 9, top: 4, filename: "map-icons.png"};
mapIconsConfig.image = new Image();
mapIconsConfig.image.onload = draw;
mapIconsConfig.image.src = mapIconsConfig.filename;

let _lastUiState = {};
function draw() {
    let map = getMap(uiState.size);
    let noisyEdges = uiState['noisy-edges'],
        noisyFills = uiState['noisy-fills'];
    
    let canvas = document.getElementById('map');
    let ctx = canvas.getContext('2d');

    let size = Math.min(canvas.parentNode.clientWidth, canvas.parentNode.clientHeight);
    if (size != uiState.canvasSize) {
        // Don't assign to width,height if the size hasn't changed because
        // it will blank out the canvas and we'd like to reuse the previous draw
        uiState.canvasSize = size;
        canvas.style.width = size + 'px';
        canvas.style.height = size + 'px';
        size = 1024;
        if (window.devicePixelRatio && window.devicePixelRatio != 1) {
            size *= window.devicePixelRatio;
        }
        canvas.width = size;
        canvas.height = size;
    }
    
    let noise = new SimplexNoise(makeRandFloat(uiState.seed));
    let persistence = Math.pow(1/2, 1 + uiState.persistence);
    let islandShapeAmplitudes = Array.from({length: 5}, (_,octave) => Math.pow(persistence, octave));
    let biomeBias = {
        north_temperature: uiState['north-temperature'],
        south_temperature: uiState['south-temperature'],
        moisture: uiState.rainfall,
    };
    let colormap = uiState.biomes? new Colormap.Discrete() : new Colormap.Smooth();
    let queue = [];
    if ((!noisyEdges || uiState.size === 'large' || uiState.size === 'huge')
        && (_lastUiState.seed !== uiState.seed
            || _lastUiState.size !== uiState.size
            || _lastUiState.canvasSize !== uiState.canvasSize)) {
        // Noisy edges are slow enough that it'd be nice to have a
        // quick approximation drawn first, but if the last time we
        // drew something was with the same essential parameters let's
        // reuse the drawing from last time
        queue.push(() => Draw.approximateIslandShape(ctx, 1000, 1000, noise, {round: 0.5, inflate: 0.4, amplitudes: islandShapeAmplitudes.slice(0, 3)}));
        // TODO: the if() test is too convoluted; rewrite that expression
    }
    Object.assign(_lastUiState, uiState);

    queue.push(
        () => map.calculate({
            noise: noise,
            drainageSeed: uiState.variant,
            riverSeed: uiState.variant,
            biomeBias: biomeBias,
            shape: {round: 0.5, inflate: 0.4, amplitudes: islandShapeAmplitudes},
        }),
        () => {
            Draw.background(ctx, colormap);
            Draw.noisyRegions(ctx, map, colormap, noisyEdges);
            // Draw the rivers early for better user experience
            Draw.rivers(ctx, map, colormap, noisyEdges, true);
        }
    );

    for (let phase = 0; phase < 16; phase++) {
        queue.push(() => Draw.noisyEdges(ctx, map, colormap, noisyEdges, phase));
    }

    // Have to draw the rivers and coastlines again because the noisy
    // edges might overwrite them, and these should take priority over
    // the other noisy edges. Otherwise it leaves little gaps that look
    // ugly when zoomed in.
    queue.push(() => Draw.rivers(ctx, map, colormap, noisyEdges, false));
    queue.push(() => Draw.coastlines(ctx, map, colormap, noisyEdges));

    if (noisyFills) {
        queue.push(() => Draw.noisyFill(ctx, 1000, 1000, makeRandInt(12345)));
    }

    if (uiState.icons) {
        queue.push(() => Draw.regionIcons(ctx, map, mapIconsConfig, makeRandInt(uiState.variant)));
    }
    
    if (uiState.lighting) {
        queue.push(() => Draw.lighting(ctx, 1000, 1000, map));
    }

    requestAnimationFrameQueue = queue.map(
        (layer, i) => () => {
            //console.time("layer "+i);
            ctx.save();
            ctx.scale(canvas.width / 1000, canvas.height / 1000);
            layer();
            ctx.restore();
            //console.timeEnd("layer "+i);
        });

    if (!requestAnimationFrameId) {
        requestAnimationFrameId = requestAnimationFrame(requestAnimationFrameHandler);
    }
}


function initUi() {
    function oninput(element) { element.addEventListener('input', getUiState); }
    function onclick(element) { element.addEventListener('click', getUiState); }
    function onchange(element) { element.addEventListener('change', getUiState); }
    document.querySelectorAll("input[type='radio']").forEach(onclick);
    document.querySelectorAll("input[type='checkbox']").forEach(onclick);
    document.querySelectorAll("input[type='number']").forEach(onchange);
    document.querySelectorAll("input[type='range']").forEach(oninput);

    // HACK: on touch devices use touch event to make the slider feel better
    document.querySelectorAll("input[type='range']").forEach((slider) => {
        function handleTouch(e) {
            let rect = slider.getBoundingClientRect();
            let min = parseFloat(slider.getAttribute('min')),
                max = parseFloat(slider.getAttribute('max')),
                step = parseFloat(slider.getAttribute('step')) || 1;
            let value = (e.changedTouches[0].clientX - rect.left) / rect.width;
            value = min + value * (max - min);
            value = Math.round(value / step) * step;
            if (value < min) { value = min; }
            if (value > max) { value = max; }
            slider.value = value;
            slider.dispatchEvent(new Event('input'));
            e.preventDefault();
            e.stopPropagation();
        };
        slider.addEventListener('touchmove', handleTouch, {passive: true});
        slider.addEventListener('touchstart', handleTouch, {passive: true});
    });
}

function setUiState() {
    document.getElementById('seed').value = uiState.seed;
    document.getElementById('variant').value = uiState.variant;
    document.querySelector("input#size-" + uiState.size).checked = true;
    document.querySelector("input#noisy-edges").checked = uiState['noisy-edges'];
    document.querySelector("input#noisy-fills").checked = uiState['noisy-fills'];
    document.querySelector("input#icons").checked = uiState.icons;
    document.querySelector("input#biomes").checked = uiState.biomes;
    document.querySelector("input#lighting").checked = uiState.lighting;
    document.querySelector("input#north-temperature").value = uiState['north-temperature'];
    document.querySelector("input#south-temperature").value = uiState['south-temperature'];
    document.querySelector("input#rainfall").value = uiState.rainfall;
    document.querySelector("input#persistence").value = uiState.persistence;
}

function getUiState() {
    uiState.seed = document.getElementById('seed').valueAsNumber;
    uiState.variant = document.getElementById('variant').valueAsNumber;
    uiState.size = document.querySelector("input[name='size']:checked").value;
    uiState['noisy-edges'] = document.querySelector("input#noisy-edges").checked;
    uiState['noisy-fills'] = document.querySelector("input#noisy-fills").checked;
    uiState.icons = document.querySelector("input#icons").checked;
    uiState.biomes = document.querySelector("input#biomes").checked;
    uiState.lighting = document.querySelector("input#lighting").checked;
    uiState['north-temperature'] = document.querySelector("input#north-temperature").valueAsNumber;
    uiState['south-temperature'] = document.querySelector("input#south-temperature").valueAsNumber;
    uiState.rainfall = document.querySelector("input#rainfall").valueAsNumber;
    uiState.persistence = document.querySelector("input#persistence").valueAsNumber;
    setUrlFromState();
    draw();
}

function setSeed(seed) {
    uiState.seed = seed & 0x7fffffff;
    setUiState();
    getUiState();
}

function setVariant(variant) {
    uiState.variant = ((variant % 10) + 10) % 10;
    setUiState();
    getUiState();
}

global.prevSeed = function() { setSeed(uiState.seed - 1); };
global.nextSeed = function() { setSeed(uiState.seed + 1); };
global.prevVariant = function() { setVariant(uiState.variant - 1); };
global.nextVariant = function() { setVariant(uiState.variant + 1); };


let _setUrlFromStateTimeout = null;
function _setUrlFromState() {
    _setUrlFromStateTimeout = null;
    let fragment = urlUtils.stringifyParams(uiState, {}, {
        'canvasSize': 'exclude',
        'noisy-edges': 'include-if-falsy',
        'noisy-fills': 'include-if-falsy',
    });
    let url = window.location.pathname + "#" + fragment;
    window.history.replaceState({}, null, url);
    document.getElementById('url').setAttribute('href', window.location.href);
}
function setUrlFromState() {
    // Rate limit the url update because some browsers (like Safari
    // iOS) throw an error if you change the url too quickly.
    if (_setUrlFromStateTimeout === null) {
        _setUrlFromStateTimeout = setTimeout(_setUrlFromState, 500);
    }
}
    
function getStateFromUrl() {
    const bool = (value) => value === 'true';
    let hashState = urlUtils.parseQuery(
        window.location.hash.slice(1),
        {
            'seed': 'number',
            'variant': 'number',
            'noisy-edges': bool,
            'noisy-fills': bool,
            'icons': bool,
            'biomes': bool,
            'lighting': bool,
            'north-temperature': 'number',
            'south-temperature': 'number',
            'rainfall': 'number',
            'persistence': 'number',
        }
    );
    Object.assign(uiState, defaultUiState);
    Object.assign(uiState, hashState);
    if (hashState.temperature) {
        // backwards compatibility -- I changed the url from
        // &temperature= to separate north and south parameters
        uiState['north-temperature'] = uiState.temperature;
        uiState['south-temperature'] = uiState.temperature;
        delete uiState.temperature;
    }
    setUrlFromState();
    setUiState();
    draw();
}
window.addEventListener('hashchange', getStateFromUrl);
window.addEventListener('resize', draw);

initUi();
getStateFromUrl();
setUiState();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./colormap":23,"./draw":24,"@redblobgames/dual-mesh":2,"@redblobgames/dual-mesh/create":1,"@redblobgames/mapgen2":15,"@redblobgames/prng":21,"simplex-noise":26,"url-search-utils":27}],26:[function(require,module,exports){
/*
 * A fast javascript implementation of simplex noise by Jonas Wagner
 *
 * Based on a speed-improved simplex noise algorithm for 2D, 3D and 4D in Java.
 * Which is based on example code by Stefan Gustavson (stegu@itn.liu.se).
 * With Optimisations by Peter Eastman (peastman@drizzle.stanford.edu).
 * Better rank ordering method by Stefan Gustavson in 2012.
 *
 *
 * Copyright (C) 2016 Jonas Wagner
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */
(function() {
'use strict';

var F2 = 0.5 * (Math.sqrt(3.0) - 1.0);
var G2 = (3.0 - Math.sqrt(3.0)) / 6.0;
var F3 = 1.0 / 3.0;
var G3 = 1.0 / 6.0;
var F4 = (Math.sqrt(5.0) - 1.0) / 4.0;
var G4 = (5.0 - Math.sqrt(5.0)) / 20.0;

function SimplexNoise(random) {
  if (!random) random = Math.random;
  this.p = buildPermutationTable(random);
  this.perm = new Uint8Array(512);
  this.permMod12 = new Uint8Array(512);
  for (var i = 0; i < 512; i++) {
    this.perm[i] = this.p[i & 255];
    this.permMod12[i] = this.perm[i] % 12;
  }

}
SimplexNoise.prototype = {
    grad3: new Float32Array([1, 1, 0,
                            -1, 1, 0,
                            1, -1, 0,

                            -1, -1, 0,
                            1, 0, 1,
                            -1, 0, 1,

                            1, 0, -1,
                            -1, 0, -1,
                            0, 1, 1,

                            0, -1, 1,
                            0, 1, -1,
                            0, -1, -1]),
    grad4: new Float32Array([0, 1, 1, 1, 0, 1, 1, -1, 0, 1, -1, 1, 0, 1, -1, -1,
                            0, -1, 1, 1, 0, -1, 1, -1, 0, -1, -1, 1, 0, -1, -1, -1,
                            1, 0, 1, 1, 1, 0, 1, -1, 1, 0, -1, 1, 1, 0, -1, -1,
                            -1, 0, 1, 1, -1, 0, 1, -1, -1, 0, -1, 1, -1, 0, -1, -1,
                            1, 1, 0, 1, 1, 1, 0, -1, 1, -1, 0, 1, 1, -1, 0, -1,
                            -1, 1, 0, 1, -1, 1, 0, -1, -1, -1, 0, 1, -1, -1, 0, -1,
                            1, 1, 1, 0, 1, 1, -1, 0, 1, -1, 1, 0, 1, -1, -1, 0,
                            -1, 1, 1, 0, -1, 1, -1, 0, -1, -1, 1, 0, -1, -1, -1, 0]),
    noise2D: function(xin, yin) {
        var permMod12 = this.permMod12;
        var perm = this.perm;
        var grad3 = this.grad3;
        var n0 = 0; // Noise contributions from the three corners
        var n1 = 0;
        var n2 = 0;
        // Skew the input space to determine which simplex cell we're in
        var s = (xin + yin) * F2; // Hairy factor for 2D
        var i = Math.floor(xin + s);
        var j = Math.floor(yin + s);
        var t = (i + j) * G2;
        var X0 = i - t; // Unskew the cell origin back to (x,y) space
        var Y0 = j - t;
        var x0 = xin - X0; // The x,y distances from the cell origin
        var y0 = yin - Y0;
        // For the 2D case, the simplex shape is an equilateral triangle.
        // Determine which simplex we are in.
        var i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
        if (x0 > y0) {
          i1 = 1;
          j1 = 0;
        } // lower triangle, XY order: (0,0)->(1,0)->(1,1)
        else {
          i1 = 0;
          j1 = 1;
        } // upper triangle, YX order: (0,0)->(0,1)->(1,1)
        // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
        // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
        // c = (3-sqrt(3))/6
        var x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
        var y1 = y0 - j1 + G2;
        var x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
        var y2 = y0 - 1.0 + 2.0 * G2;
        // Work out the hashed gradient indices of the three simplex corners
        var ii = i & 255;
        var jj = j & 255;
        // Calculate the contribution from the three corners
        var t0 = 0.5 - x0 * x0 - y0 * y0;
        if (t0 >= 0) {
          var gi0 = permMod12[ii + perm[jj]] * 3;
          t0 *= t0;
          n0 = t0 * t0 * (grad3[gi0] * x0 + grad3[gi0 + 1] * y0); // (x,y) of grad3 used for 2D gradient
        }
        var t1 = 0.5 - x1 * x1 - y1 * y1;
        if (t1 >= 0) {
          var gi1 = permMod12[ii + i1 + perm[jj + j1]] * 3;
          t1 *= t1;
          n1 = t1 * t1 * (grad3[gi1] * x1 + grad3[gi1 + 1] * y1);
        }
        var t2 = 0.5 - x2 * x2 - y2 * y2;
        if (t2 >= 0) {
          var gi2 = permMod12[ii + 1 + perm[jj + 1]] * 3;
          t2 *= t2;
          n2 = t2 * t2 * (grad3[gi2] * x2 + grad3[gi2 + 1] * y2);
        }
        // Add contributions from each corner to get the final noise value.
        // The result is scaled to return values in the interval [-1,1].
        return 70.0 * (n0 + n1 + n2);
      },
    // 3D simplex noise
    noise3D: function(xin, yin, zin) {
        var permMod12 = this.permMod12;
        var perm = this.perm;
        var grad3 = this.grad3;
        var n0, n1, n2, n3; // Noise contributions from the four corners
        // Skew the input space to determine which simplex cell we're in
        var s = (xin + yin + zin) * F3; // Very nice and simple skew factor for 3D
        var i = Math.floor(xin + s);
        var j = Math.floor(yin + s);
        var k = Math.floor(zin + s);
        var t = (i + j + k) * G3;
        var X0 = i - t; // Unskew the cell origin back to (x,y,z) space
        var Y0 = j - t;
        var Z0 = k - t;
        var x0 = xin - X0; // The x,y,z distances from the cell origin
        var y0 = yin - Y0;
        var z0 = zin - Z0;
        // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
        // Determine which simplex we are in.
        var i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
        var i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
        if (x0 >= y0) {
          if (y0 >= z0) {
            i1 = 1;
            j1 = 0;
            k1 = 0;
            i2 = 1;
            j2 = 1;
            k2 = 0;
          } // X Y Z order
          else if (x0 >= z0) {
            i1 = 1;
            j1 = 0;
            k1 = 0;
            i2 = 1;
            j2 = 0;
            k2 = 1;
          } // X Z Y order
          else {
            i1 = 0;
            j1 = 0;
            k1 = 1;
            i2 = 1;
            j2 = 0;
            k2 = 1;
          } // Z X Y order
        }
        else { // x0<y0
          if (y0 < z0) {
            i1 = 0;
            j1 = 0;
            k1 = 1;
            i2 = 0;
            j2 = 1;
            k2 = 1;
          } // Z Y X order
          else if (x0 < z0) {
            i1 = 0;
            j1 = 1;
            k1 = 0;
            i2 = 0;
            j2 = 1;
            k2 = 1;
          } // Y Z X order
          else {
            i1 = 0;
            j1 = 1;
            k1 = 0;
            i2 = 1;
            j2 = 1;
            k2 = 0;
          } // Y X Z order
        }
        // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
        // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
        // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
        // c = 1/6.
        var x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
        var y1 = y0 - j1 + G3;
        var z1 = z0 - k1 + G3;
        var x2 = x0 - i2 + 2.0 * G3; // Offsets for third corner in (x,y,z) coords
        var y2 = y0 - j2 + 2.0 * G3;
        var z2 = z0 - k2 + 2.0 * G3;
        var x3 = x0 - 1.0 + 3.0 * G3; // Offsets for last corner in (x,y,z) coords
        var y3 = y0 - 1.0 + 3.0 * G3;
        var z3 = z0 - 1.0 + 3.0 * G3;
        // Work out the hashed gradient indices of the four simplex corners
        var ii = i & 255;
        var jj = j & 255;
        var kk = k & 255;
        // Calculate the contribution from the four corners
        var t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
        if (t0 < 0) n0 = 0.0;
        else {
          var gi0 = permMod12[ii + perm[jj + perm[kk]]] * 3;
          t0 *= t0;
          n0 = t0 * t0 * (grad3[gi0] * x0 + grad3[gi0 + 1] * y0 + grad3[gi0 + 2] * z0);
        }
        var t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
        if (t1 < 0) n1 = 0.0;
        else {
          var gi1 = permMod12[ii + i1 + perm[jj + j1 + perm[kk + k1]]] * 3;
          t1 *= t1;
          n1 = t1 * t1 * (grad3[gi1] * x1 + grad3[gi1 + 1] * y1 + grad3[gi1 + 2] * z1);
        }
        var t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
        if (t2 < 0) n2 = 0.0;
        else {
          var gi2 = permMod12[ii + i2 + perm[jj + j2 + perm[kk + k2]]] * 3;
          t2 *= t2;
          n2 = t2 * t2 * (grad3[gi2] * x2 + grad3[gi2 + 1] * y2 + grad3[gi2 + 2] * z2);
        }
        var t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
        if (t3 < 0) n3 = 0.0;
        else {
          var gi3 = permMod12[ii + 1 + perm[jj + 1 + perm[kk + 1]]] * 3;
          t3 *= t3;
          n3 = t3 * t3 * (grad3[gi3] * x3 + grad3[gi3 + 1] * y3 + grad3[gi3 + 2] * z3);
        }
        // Add contributions from each corner to get the final noise value.
        // The result is scaled to stay just inside [-1,1]
        return 32.0 * (n0 + n1 + n2 + n3);
      },
    // 4D simplex noise, better simplex rank ordering method 2012-03-09
    noise4D: function(x, y, z, w) {
        var permMod12 = this.permMod12;
        var perm = this.perm;
        var grad4 = this.grad4;

        var n0, n1, n2, n3, n4; // Noise contributions from the five corners
        // Skew the (x,y,z,w) space to determine which cell of 24 simplices we're in
        var s = (x + y + z + w) * F4; // Factor for 4D skewing
        var i = Math.floor(x + s);
        var j = Math.floor(y + s);
        var k = Math.floor(z + s);
        var l = Math.floor(w + s);
        var t = (i + j + k + l) * G4; // Factor for 4D unskewing
        var X0 = i - t; // Unskew the cell origin back to (x,y,z,w) space
        var Y0 = j - t;
        var Z0 = k - t;
        var W0 = l - t;
        var x0 = x - X0; // The x,y,z,w distances from the cell origin
        var y0 = y - Y0;
        var z0 = z - Z0;
        var w0 = w - W0;
        // For the 4D case, the simplex is a 4D shape I won't even try to describe.
        // To find out which of the 24 possible simplices we're in, we need to
        // determine the magnitude ordering of x0, y0, z0 and w0.
        // Six pair-wise comparisons are performed between each possible pair
        // of the four coordinates, and the results are used to rank the numbers.
        var rankx = 0;
        var ranky = 0;
        var rankz = 0;
        var rankw = 0;
        if (x0 > y0) rankx++;
        else ranky++;
        if (x0 > z0) rankx++;
        else rankz++;
        if (x0 > w0) rankx++;
        else rankw++;
        if (y0 > z0) ranky++;
        else rankz++;
        if (y0 > w0) ranky++;
        else rankw++;
        if (z0 > w0) rankz++;
        else rankw++;
        var i1, j1, k1, l1; // The integer offsets for the second simplex corner
        var i2, j2, k2, l2; // The integer offsets for the third simplex corner
        var i3, j3, k3, l3; // The integer offsets for the fourth simplex corner
        // simplex[c] is a 4-vector with the numbers 0, 1, 2 and 3 in some order.
        // Many values of c will never occur, since e.g. x>y>z>w makes x<z, y<w and x<w
        // impossible. Only the 24 indices which have non-zero entries make any sense.
        // We use a thresholding to set the coordinates in turn from the largest magnitude.
        // Rank 3 denotes the largest coordinate.
        i1 = rankx >= 3 ? 1 : 0;
        j1 = ranky >= 3 ? 1 : 0;
        k1 = rankz >= 3 ? 1 : 0;
        l1 = rankw >= 3 ? 1 : 0;
        // Rank 2 denotes the second largest coordinate.
        i2 = rankx >= 2 ? 1 : 0;
        j2 = ranky >= 2 ? 1 : 0;
        k2 = rankz >= 2 ? 1 : 0;
        l2 = rankw >= 2 ? 1 : 0;
        // Rank 1 denotes the second smallest coordinate.
        i3 = rankx >= 1 ? 1 : 0;
        j3 = ranky >= 1 ? 1 : 0;
        k3 = rankz >= 1 ? 1 : 0;
        l3 = rankw >= 1 ? 1 : 0;
        // The fifth corner has all coordinate offsets = 1, so no need to compute that.
        var x1 = x0 - i1 + G4; // Offsets for second corner in (x,y,z,w) coords
        var y1 = y0 - j1 + G4;
        var z1 = z0 - k1 + G4;
        var w1 = w0 - l1 + G4;
        var x2 = x0 - i2 + 2.0 * G4; // Offsets for third corner in (x,y,z,w) coords
        var y2 = y0 - j2 + 2.0 * G4;
        var z2 = z0 - k2 + 2.0 * G4;
        var w2 = w0 - l2 + 2.0 * G4;
        var x3 = x0 - i3 + 3.0 * G4; // Offsets for fourth corner in (x,y,z,w) coords
        var y3 = y0 - j3 + 3.0 * G4;
        var z3 = z0 - k3 + 3.0 * G4;
        var w3 = w0 - l3 + 3.0 * G4;
        var x4 = x0 - 1.0 + 4.0 * G4; // Offsets for last corner in (x,y,z,w) coords
        var y4 = y0 - 1.0 + 4.0 * G4;
        var z4 = z0 - 1.0 + 4.0 * G4;
        var w4 = w0 - 1.0 + 4.0 * G4;
        // Work out the hashed gradient indices of the five simplex corners
        var ii = i & 255;
        var jj = j & 255;
        var kk = k & 255;
        var ll = l & 255;
        // Calculate the contribution from the five corners
        var t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0 - w0 * w0;
        if (t0 < 0) n0 = 0.0;
        else {
          var gi0 = (perm[ii + perm[jj + perm[kk + perm[ll]]]] % 32) * 4;
          t0 *= t0;
          n0 = t0 * t0 * (grad4[gi0] * x0 + grad4[gi0 + 1] * y0 + grad4[gi0 + 2] * z0 + grad4[gi0 + 3] * w0);
        }
        var t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1 - w1 * w1;
        if (t1 < 0) n1 = 0.0;
        else {
          var gi1 = (perm[ii + i1 + perm[jj + j1 + perm[kk + k1 + perm[ll + l1]]]] % 32) * 4;
          t1 *= t1;
          n1 = t1 * t1 * (grad4[gi1] * x1 + grad4[gi1 + 1] * y1 + grad4[gi1 + 2] * z1 + grad4[gi1 + 3] * w1);
        }
        var t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2 - w2 * w2;
        if (t2 < 0) n2 = 0.0;
        else {
          var gi2 = (perm[ii + i2 + perm[jj + j2 + perm[kk + k2 + perm[ll + l2]]]] % 32) * 4;
          t2 *= t2;
          n2 = t2 * t2 * (grad4[gi2] * x2 + grad4[gi2 + 1] * y2 + grad4[gi2 + 2] * z2 + grad4[gi2 + 3] * w2);
        }
        var t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3 - w3 * w3;
        if (t3 < 0) n3 = 0.0;
        else {
          var gi3 = (perm[ii + i3 + perm[jj + j3 + perm[kk + k3 + perm[ll + l3]]]] % 32) * 4;
          t3 *= t3;
          n3 = t3 * t3 * (grad4[gi3] * x3 + grad4[gi3 + 1] * y3 + grad4[gi3 + 2] * z3 + grad4[gi3 + 3] * w3);
        }
        var t4 = 0.6 - x4 * x4 - y4 * y4 - z4 * z4 - w4 * w4;
        if (t4 < 0) n4 = 0.0;
        else {
          var gi4 = (perm[ii + 1 + perm[jj + 1 + perm[kk + 1 + perm[ll + 1]]]] % 32) * 4;
          t4 *= t4;
          n4 = t4 * t4 * (grad4[gi4] * x4 + grad4[gi4 + 1] * y4 + grad4[gi4 + 2] * z4 + grad4[gi4 + 3] * w4);
        }
        // Sum up and scale the result to cover the range [-1,1]
        return 27.0 * (n0 + n1 + n2 + n3 + n4);
      }
  };

function buildPermutationTable(random) {
  var i;
  var p = new Uint8Array(256);
  for (i = 0; i < 256; i++) {
    p[i] = i;
  }
  for (i = 0; i < 255; i++) {
    var r = i + 1 + ~~(random() * (255 - i));
    var aux = p[i];
    p[i] = p[r];
    p[r] = aux;
  }
  return p;
}
SimplexNoise._buildPermutationTable = buildPermutationTable;

// amd
if (typeof define !== 'undefined' && define.amd) define(function() {return SimplexNoise;});
// common js
if (typeof exports !== 'undefined') exports.SimplexNoise = SimplexNoise;
// browser
else if (typeof window !== 'undefined') window.SimplexNoise = SimplexNoise;
// nodejs
if (typeof module !== 'undefined') {
  module.exports = SimplexNoise;
}

})();

},{}],27:[function(require,module,exports){
(function webpackUniversalModuleDefinition(root, factory) {
  if(typeof exports === 'object' && typeof module === 'object')
    module.exports = factory();
  else if(typeof define === 'function' && define.amd)
    define([], factory);
  else if(typeof exports === 'object')
    exports["URLSearchUtils"] = factory();
  else
    root["URLSearchUtils"] = factory();
})(this, function() {
return /******/ (function(modules) { // webpackBootstrap
/******/  // The module cache
/******/  var installedModules = {};
/******/
/******/  // The require function
/******/  function __webpack_require__(moduleId) {
/******/
/******/    // Check if module is in cache
/******/    if(installedModules[moduleId]) {
/******/      return installedModules[moduleId].exports;
/******/    }
/******/    // Create a new module (and put it into the cache)
/******/    var module = installedModules[moduleId] = {
/******/      i: moduleId,
/******/      l: false,
/******/      exports: {}
/******/    };
/******/
/******/    // Execute the module function
/******/    modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/
/******/    // Flag the module as loaded
/******/    module.l = true;
/******/
/******/    // Return the exports of the module
/******/    return module.exports;
/******/  }
/******/
/******/
/******/  // expose the modules object (__webpack_modules__)
/******/  __webpack_require__.m = modules;
/******/
/******/  // expose the module cache
/******/  __webpack_require__.c = installedModules;
/******/
/******/  // identity function for calling harmony imports with the correct context
/******/  __webpack_require__.i = function(value) { return value; };
/******/
/******/  // define getter function for harmony exports
/******/  __webpack_require__.d = function(exports, name, getter) {
/******/    if(!__webpack_require__.o(exports, name)) {
/******/      Object.defineProperty(exports, name, {
/******/        configurable: false,
/******/        enumerable: true,
/******/        get: getter
/******/      });
/******/    }
/******/  };
/******/
/******/  // getDefaultExport function for compatibility with non-harmony modules
/******/  __webpack_require__.n = function(module) {
/******/    var getter = module && module.__esModule ?
/******/      function getDefault() { return module['default']; } :
/******/      function getModuleExports() { return module; };
/******/    __webpack_require__.d(getter, 'a', getter);
/******/    return getter;
/******/  };
/******/
/******/  // Object.prototype.hasOwnProperty.call
/******/  __webpack_require__.o = function(object, property) { return Object.prototype.hasOwnProperty.call(object, property); };
/******/
/******/  // __webpack_public_path__
/******/  __webpack_require__.p = "";
/******/
/******/  // Load entry module and return exports
/******/  return __webpack_require__(__webpack_require__.s = 6);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

exports.default = parseQuery;
function parseQuery(search) {
  var paramsTypes = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {};

  if (search.length === 0) {
    return {};
  }

  return search.split('&').reduce(function (result, searchItem) {
    var _searchItem$split = searchItem.split('='),
        _searchItem$split2 = _slicedToArray(_searchItem$split, 2),
        name = _searchItem$split2[0],
        rawValue = _searchItem$split2[1];

    var value = decodeURIComponent(rawValue);

    if (paramsTypes[name]) {
      if (typeof paramsTypes[name] === 'function') {
        result[name] = paramsTypes[name](value, result[name]);

        return result;
      }

      switch (paramsTypes[name]) {
        case 'array-of-strings':
          if (!result[name]) {
            result[name] = [];
          }

          result[name].push(value);

          return result;

        case 'array-of-numbers':
          if (!result[name]) {
            result[name] = [];
          }

          result[name].push(parseFloat(value));

          return result;

        case 'number':
          result[name] = parseFloat(value);

          return result;

        case 'exclude':

          return result;

        default:
          throw new Error('Unknown type of parameter for parse "' + name + '": "' + paramsTypes[name] + '"');
      }
    }

    result[name] = value;

    return result;
  }, {});
}

/***/ }),
/* 1 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = stringifyParams;
function stringifyParams(params) {
  var mapParamsNames = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {};
  var paramsTypes = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : {};

  return Object.keys(params).map(function (rawParamName) {
    var paramName = mapParamsNames[rawParamName] || rawParamName;
    var paramValue = params[rawParamName];

    if (paramsTypes[paramName]) {
      if (typeof paramsTypes[paramName] === 'function') {
        var mappedParamValue = paramsTypes[paramName](paramValue);

        if (mappedParamValue === null || typeof mappedParamValue === 'undefined') {
          return null;
        }

        return paramName + '=' + encodeURIComponent(mappedParamValue);
      }

      switch (paramsTypes[paramName]) {
        case 'exclude':
          return null;

        case 'include-if-falsy':
          return paramName + '=' + encodeURIComponent(paramValue || '');

        default:
          throw new Error('Unknown type of parameter for serialize "' + paramName + '": "' + paramsTypes[paramName] + '"');
      }
    }

    if (paramValue instanceof Array) {
      if (paramValue.length === 0) {
        return null;
      }

      return paramValue.map(function (paramValueItem) {
        return paramName + '=' + encodeURIComponent(paramValueItem);
      }).join('&');
    }

    if (!paramValue) {
      return null;
    }

    return paramName + '=' + encodeURIComponent(paramValue);
  }).filter(function (item) {
    return item !== null;
  }).join('&');
}

/***/ }),
/* 2 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = parseHashParams;

var _parseQuery = __webpack_require__(0);

var _parseQuery2 = _interopRequireDefault(_parseQuery);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function parseHashParams() {
  var paramsTypes = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};

  var search = window.location.hash.split('?')[1] || '';

  return (0, _parseQuery2.default)(search, paramsTypes);
}

/***/ }),
/* 3 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = parseSearchParams;

var _parseQuery = __webpack_require__(0);

var _parseQuery2 = _interopRequireDefault(_parseQuery);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function parseSearchParams() {
  var paramsTypes = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};

  var search = window.location.search;

  return (0, _parseQuery2.default)(search.substring(1, search.length), paramsTypes);
}

/***/ }),
/* 4 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = setHashParams;

var _stringifyParams = __webpack_require__(1);

var _stringifyParams2 = _interopRequireDefault(_stringifyParams);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function setHashParams(params, mapParamsNames, paramsTypes) {
  var paramsStr = (0, _stringifyParams2.default)(params, mapParamsNames, paramsTypes);

  var hashWithoutSearch = window.location.hash.split('?')[0];

  window.location.hash = hashWithoutSearch + '?' + paramsStr;
}

/***/ }),
/* 5 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = setSearchParams;

var _stringifyParams = __webpack_require__(1);

var _stringifyParams2 = _interopRequireDefault(_stringifyParams);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function setSearchParams(params, mapParamsNames, paramsTypes) {
  var paramsStr = (0, _stringifyParams2.default)(params, mapParamsNames, paramsTypes);

  window.history.replaceState({}, null, window.location.pathname + '?' + paramsStr);
}

/***/ }),
/* 6 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _parseQuery = __webpack_require__(0);

Object.defineProperty(exports, 'parseQuery', {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_parseQuery).default;
  }
});

var _stringifyParams = __webpack_require__(1);

Object.defineProperty(exports, 'stringifyParams', {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_stringifyParams).default;
  }
});

var _parseSearchParams = __webpack_require__(3);

Object.defineProperty(exports, 'parseSearchParams', {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_parseSearchParams).default;
  }
});

var _parseHashParams = __webpack_require__(2);

Object.defineProperty(exports, 'parseHashParams', {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_parseHashParams).default;
  }
});

var _setSearchParams = __webpack_require__(5);

Object.defineProperty(exports, 'setSearchParams', {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_setSearchParams).default;
  }
});

var _setHashParams = __webpack_require__(4);

Object.defineProperty(exports, 'setHashParams', {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_setHashParams).default;
  }
});

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

/***/ })
/******/ ]);
});
},{}]},{},[25]);
document.write(`      <!-- Start of counters -->
      <script>
var sc_project=417499;
var sc_invisible=1;
var sc_security="";
(function() {
var script = document.createElement('script');
script.type = 'text/javascript'; script.async = true;
script.src = 'http://statcounter.com/counter/counter_xhtml.js';
if (document.location.hostname != 'localhost') document.getElementsByTagName('body')[0].appendChild(script);
})();
var _gap = _gap || [];
_gap.push(['_setAccount', 'UA-79181-1']);
_gap.push(['_setDomainName', 'redblobgames.com']);
_gap.push(['_setAllowLinker', true]);
_gap.push(['_trackPageview']);
_gap.push(['_gapTrackBounceViaTime', 30]);
_gap.push(['_gapTrackBounceViaScroll', 25]);
_gap.push(['_gapTrackReads', 60, 10]);
_gap.push(['_gapTrackLinkClicks']);
if (document.location.hostname != 'localhost') (function() {
var gap = document.createElement('script');
gap.async = true;
gap.type = 'text/javascript';
gap.src = '/js/gap.min.js';
var s = document.getElementsByTagName('script')[0];
s.parentNode.insertBefore(gap, s);
})();
       </script>
       <!-- End of counters -->      `);
       
