
// const EPSILON = Math.pow(2, -52);
// const EDGE_STACK = new Uint32Array(512);

//  class Delaunator {

//     static from(points, getX = defaultGetX, getY = defaultGetY) {
//         const n = points.length;
//         const coords = new Float64Array(n * 2);

//         for (let i = 0; i < n; i++) {
//             const p = points[i];
//             coords[2 * i] = getX(p);
//             coords[2 * i + 1] = getY(p);
//         }

//         return new Delaunator(coords);
//     }

//     constructor(coords) {
//         const n = coords.length >> 1;
//         if (n > 0 && typeof coords[0] !== 'number') throw new Error('Expected coords to contain numbers.');

//         this.coords = coords;

//         // arrays that will store the triangulation graph
//         const maxTriangles = Math.max(2 * n - 5, 0);
//         this._triangles = new Uint32Array(maxTriangles * 3);
//         this._halfedges = new Int32Array(maxTriangles * 3);

//         // temporary arrays for tracking the edges of the advancing convex hull
//         this._hashSize = Math.ceil(Math.sqrt(n));
//         this._hullPrev = new Uint32Array(n); // edge to prev edge
//         this._hullNext = new Uint32Array(n); // edge to next edge
//         this._hullTri = new Uint32Array(n); // edge to adjacent triangle
//         this._hullHash = new Int32Array(this._hashSize).fill(-1); // angular edge hash

//         // temporary arrays for sorting points
//         this._ids = new Uint32Array(n);
//         this._dists = new Float64Array(n);

//         this.update();
//     }

//     update() {
//         const {coords, _hullPrev: hullPrev, _hullNext: hullNext, _hullTri: hullTri, _hullHash: hullHash} =  this;
//         const n = coords.length >> 1;

//         // populate an array of point indices; calculate input data bbox
//         let minX = Infinity;
//         let minY = Infinity;
//         let maxX = -Infinity;
//         let maxY = -Infinity;

//         for (let i = 0; i < n; i++) {
//             const x = coords[2 * i];
//             const y = coords[2 * i + 1];
//             if (x < minX) minX = x;
//             if (y < minY) minY = y;
//             if (x > maxX) maxX = x;
//             if (y > maxY) maxY = y;
//             this._ids[i] = i;
//         }
//         const cx = (minX + maxX) / 2;
//         const cy = (minY + maxY) / 2;

//         let minDist = Infinity;
//         let i0, i1, i2;

//         // pick a seed point close to the center
//         for (let i = 0; i < n; i++) {
//             const d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
//             if (d < minDist) {
//                 i0 = i;
//                 minDist = d;
//             }
//         }
//         const i0x = coords[2 * i0];
//         const i0y = coords[2 * i0 + 1];

//         minDist = Infinity;

//         // find the point closest to the seed
//         for (let i = 0; i < n; i++) {
//             if (i === i0) continue;
//             const d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
//             if (d < minDist && d > 0) {
//                 i1 = i;
//                 minDist = d;
//             }
//         }
//         let i1x = coords[2 * i1];
//         let i1y = coords[2 * i1 + 1];

//         let minRadius = Infinity;

//         // find the third point which forms the smallest circumcircle with the first two
//         for (let i = 0; i < n; i++) {
//             if (i === i0 || i === i1) continue;
//             const r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);
//             if (r < minRadius) {
//                 i2 = i;
//                 minRadius = r;
//             }
//         }
//         let i2x = coords[2 * i2];
//         let i2y = coords[2 * i2 + 1];

//         if (minRadius === Infinity) {
//             // order collinear points by dx (or dy if all x are identical)
//             // and return the list as a hull
//             for (let i = 0; i < n; i++) {
//                 this._dists[i] = (coords[2 * i] - coords[0]) || (coords[2 * i + 1] - coords[1]);
//             }
//             quicksort(this._ids, this._dists, 0, n - 1);
//             const hull = new Uint32Array(n);
//             let j = 0;
//             for (let i = 0, d0 = -Infinity; i < n; i++) {
//                 const id = this._ids[i];
//                 if (this._dists[id] > d0) {
//                     hull[j++] = id;
//                     d0 = this._dists[id];
//                 }
//             }
//             this.hull = hull.subarray(0, j);
//             this.triangles = new Uint32Array(0);
//             this.halfedges = new Uint32Array(0);
//             return;
//         }

//         // swap the order of the seed points for counter-clockwise orientation
//         if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
//             const i = i1;
//             const x = i1x;
//             const y = i1y;
//             i1 = i2;
//             i1x = i2x;
//             i1y = i2y;
//             i2 = i;
//             i2x = x;
//             i2y = y;
//         }

//         const center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
//         this._cx = center.x;
//         this._cy = center.y;

//         for (let i = 0; i < n; i++) {
//             this._dists[i] = dist(coords[2 * i], coords[2 * i + 1], center.x, center.y);
//         }

//         // sort the points by distance from the seed triangle circumcenter
//         quicksort(this._ids, this._dists, 0, n - 1);

//         // set up the seed triangle as the starting hull
//         this._hullStart = i0;
//         let hullSize = 3;

//         hullNext[i0] = hullPrev[i2] = i1;
//         hullNext[i1] = hullPrev[i0] = i2;
//         hullNext[i2] = hullPrev[i1] = i0;

//         hullTri[i0] = 0;
//         hullTri[i1] = 1;
//         hullTri[i2] = 2;

//         hullHash.fill(-1);
//         hullHash[this._hashKey(i0x, i0y)] = i0;
//         hullHash[this._hashKey(i1x, i1y)] = i1;
//         hullHash[this._hashKey(i2x, i2y)] = i2;

//         this.trianglesLen = 0;
//         this._addTriangle(i0, i1, i2, -1, -1, -1);

//         for (let k = 0, xp, yp; k < this._ids.length; k++) {
//             const i = this._ids[k];
//             const x = coords[2 * i];
//             const y = coords[2 * i + 1];

//             // skip near-duplicate points
//             if (k > 0 && Math.abs(x - xp) <= EPSILON && Math.abs(y - yp) <= EPSILON) continue;
//             xp = x;
//             yp = y;

//             // skip seed triangle points
//             if (i === i0 || i === i1 || i === i2) continue;

//             // find a visible edge on the convex hull using edge hash
//             let start = 0;
//             for (let j = 0, key = this._hashKey(x, y); j < this._hashSize; j++) {
//                 start = hullHash[(key + j) % this._hashSize];
//                 if (start !== -1 && start !== hullNext[start]) break;
//             }

//             start = hullPrev[start];
//             let e = start, q;
//             while (q = hullNext[e], !orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
//                 e = q;
//                 if (e === start) {
//                     e = -1;
//                     break;
//                 }
//             }
//             if (e === -1) continue; // likely a near-duplicate point; skip it

//             // add the first triangle from the point
//             let t = this._addTriangle(e, i, hullNext[e], -1, -1, hullTri[e]);

//             // recursively flip triangles from the point until they satisfy the Delaunay condition
//             hullTri[i] = this._legalize(t + 2);
//             hullTri[e] = t; // keep track of boundary triangles on the hull
//             hullSize++;

//             // walk forward through the hull, adding more triangles and flipping recursively
//             let n = hullNext[e];
//             while (q = hullNext[n], orient(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1])) {
//                 t = this._addTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
//                 hullTri[i] = this._legalize(t + 2);
//                 hullNext[n] = n; // mark as removed
//                 hullSize--;
//                 n = q;
//             }

//             // walk backward from the other side, adding more triangles and flipping
//             if (e === start) {
//                 while (q = hullPrev[e], orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
//                     t = this._addTriangle(q, i, e, -1, hullTri[e], hullTri[q]);
//                     this._legalize(t + 2);
//                     hullTri[q] = t;
//                     hullNext[e] = e; // mark as removed
//                     hullSize--;
//                     e = q;
//                 }
//             }

//             // update the hull indices
//             this._hullStart = hullPrev[i] = e;
//             hullNext[e] = hullPrev[n] = i;
//             hullNext[i] = n;

//             // save the two new edges in the hash table
//             hullHash[this._hashKey(x, y)] = i;
//             hullHash[this._hashKey(coords[2 * e], coords[2 * e + 1])] = e;
//         }

//         this.hull = new Uint32Array(hullSize);
//         for (let i = 0, e = this._hullStart; i < hullSize; i++) {
//             this.hull[i] = e;
//             e = hullNext[e];
//         }

//         // trim typed triangle mesh arrays
//         this.triangles = this._triangles.subarray(0, this.trianglesLen);
//         this.halfedges = this._halfedges.subarray(0, this.trianglesLen);
//     }

//     _hashKey(x, y) {
//         return Math.floor(pseudoAngle(x - this._cx, y - this._cy) * this._hashSize) % this._hashSize;
//     }

//     _legalize(a) {
//         const {_triangles: triangles, _halfedges: halfedges, coords} = this;

//         let i = 0;
//         let ar = 0;

//         // recursion eliminated with a fixed-size stack
//         while (true) {
//             const b = halfedges[a];

//             /* if the pair of triangles doesn't satisfy the Delaunay condition
//              * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
//              * then do the same check/flip recursively for the new pair of triangles
//              *
//              *           pl                    pl
//              *          /||\                  /  \
//              *       al/ || \bl            al/    \a
//              *        /  ||  \              /      \
//              *       /  a||b  \    flip    /___ar___\
//              *     p0\   ||   /p1   =>   p0\---bl---/p1
//              *        \  ||  /              \      /
//              *       ar\ || /br             b\    /br
//              *          \||/                  \  /
//              *           pr                    pr
//              */
//             const a0 = a - a % 3;
//             ar = a0 + (a + 2) % 3;

//             if (b === -1) { // convex hull edge
//                 if (i === 0) break;
//                 a = EDGE_STACK[--i];
//                 continue;
//             }

//             const b0 = b - b % 3;
//             const al = a0 + (a + 1) % 3;
//             const bl = b0 + (b + 2) % 3;

//             const p0 = triangles[ar];
//             const pr = triangles[a];
//             const pl = triangles[al];
//             const p1 = triangles[bl];

//             const illegal = inCircle(
//                 coords[2 * p0], coords[2 * p0 + 1],
//                 coords[2 * pr], coords[2 * pr + 1],
//                 coords[2 * pl], coords[2 * pl + 1],
//                 coords[2 * p1], coords[2 * p1 + 1]);

//             if (illegal) {
//                 triangles[a] = p1;
//                 triangles[b] = p0;

//                 const hbl = halfedges[bl];

//                 // edge swapped on the other side of the hull (rare); fix the halfedge reference
//                 if (hbl === -1) {
//                     let e = this._hullStart;
//                     do {
//                         if (this._hullTri[e] === bl) {
//                             this._hullTri[e] = a;
//                             break;
//                         }
//                         e = this._hullPrev[e];
//                     } while (e !== this._hullStart);
//                 }
//                 this._link(a, hbl);
//                 this._link(b, halfedges[ar]);
//                 this._link(ar, bl);

//                 const br = b0 + (b + 1) % 3;

//                 // don't worry about hitting the cap: it can only happen on extremely degenerate input
//                 if (i < EDGE_STACK.length) {
//                     EDGE_STACK[i++] = br;
//                 }
//             } else {
//                 if (i === 0) break;
//                 a = EDGE_STACK[--i];
//             }
//         }

//         return ar;
//     }

//     _link(a, b) {
//         this._halfedges[a] = b;
//         if (b !== -1) this._halfedges[b] = a;
//     }

//     // add a new triangle given vertex indices and adjacent half-edge ids
//     _addTriangle(i0, i1, i2, a, b, c) {
//         const t = this.trianglesLen;

//         this._triangles[t] = i0;
//         this._triangles[t + 1] = i1;
//         this._triangles[t + 2] = i2;

//         this._link(t, a);
//         this._link(t + 1, b);
//         this._link(t + 2, c);

//         this.trianglesLen += 3;

//         return t;
//     }
// }

// // monotonically increases with real angle, but doesn't need expensive trigonometry
// function pseudoAngle(dx, dy) {
//     const p = dx / (Math.abs(dx) + Math.abs(dy));
//     return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
// }

// function dist(ax, ay, bx, by) {
//     const dx = ax - bx;
//     const dy = ay - by;
//     return dx * dx + dy * dy;
// }

// // return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
// function orientIfSure(px, py, rx, ry, qx, qy) {
//     const l = (ry - py) * (qx - px);
//     const r = (rx - px) * (qy - py);
//     return Math.abs(l - r) >= 3.3306690738754716e-16 * Math.abs(l + r) ? l - r : 0;
// }

// // a more robust orientation test that's stable in a given triangle (to fix robustness issues)
// function orient(rx, ry, qx, qy, px, py) {
//     return (orientIfSure(px, py, rx, ry, qx, qy) ||
//         orientIfSure(rx, ry, qx, qy, px, py) ||
//         orientIfSure(qx, qy, px, py, rx, ry)) < 0;
// }

// function inCircle(ax, ay, bx, by, cx, cy, px, py) {
//     const dx = ax - px;
//     const dy = ay - py;
//     const ex = bx - px;
//     const ey = by - py;
//     const fx = cx - px;
//     const fy = cy - py;

//     const ap = dx * dx + dy * dy;
//     const bp = ex * ex + ey * ey;
//     const cp = fx * fx + fy * fy;

//     return dx * (ey * cp - bp * fy) -
//            dy * (ex * cp - bp * fx) +
//            ap * (ex * fy - ey * fx) < 0;
// }

// function circumradius(ax, ay, bx, by, cx, cy) {
//     const dx = bx - ax;
//     const dy = by - ay;
//     const ex = cx - ax;
//     const ey = cy - ay;

//     const bl = dx * dx + dy * dy;
//     const cl = ex * ex + ey * ey;
//     const d = 0.5 / (dx * ey - dy * ex);

//     const x = (ey * bl - dy * cl) * d;
//     const y = (dx * cl - ex * bl) * d;

//     return x * x + y * y;
// }

// function circumcenter(ax, ay, bx, by, cx, cy) {
//     const dx = bx - ax;
//     const dy = by - ay;
//     const ex = cx - ax;
//     const ey = cy - ay;

//     const bl = dx * dx + dy * dy;
//     const cl = ex * ex + ey * ey;
//     const d = 0.5 / (dx * ey - dy * ex);

//     const x = ax + (ey * bl - dy * cl) * d;
//     const y = ay + (dx * cl - ex * bl) * d;

//     return {x, y};
// }

// function quicksort(ids, dists, left, right) {
//     if (right - left <= 20) {
//         for (let i = left + 1; i <= right; i++) {
//             const temp = ids[i];
//             const tempDist = dists[temp];
//             let j = i - 1;
//             while (j >= left && dists[ids[j]] > tempDist) ids[j + 1] = ids[j--];
//             ids[j + 1] = temp;
//         }
//     } else {
//         const median = (left + right) >> 1;
//         let i = left + 1;
//         let j = right;
//         swap(ids, median, i);
//         if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right);
//         if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right);
//         if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i);

//         const temp = ids[i];
//         const tempDist = dists[temp];
//         while (true) {
//             do i++; while (dists[ids[i]] < tempDist);
//             do j--; while (dists[ids[j]] > tempDist);
//             if (j < i) break;
//             swap(ids, i, j);
//         }
//         ids[left + 1] = ids[j];
//         ids[j] = temp;

//         if (right - i + 1 >= j - left) {
//             quicksort(ids, dists, i, right);
//             quicksort(ids, dists, left, j - 1);
//         } else {
//             quicksort(ids, dists, left, j - 1);
//             quicksort(ids, dists, i, right);
//         }
//     }
// }

// function swap(arr, i, j) {
//     const tmp = arr[i];
//     arr[i] = arr[j];
//     arr[j] = tmp;
// }

// function defaultGetX(p) {
//     return p[0];
// }
// function defaultGetY(p) {
//     return p[1];
// }

//eww, var
var Delaunay;

(function() {
  "use strict";

  var EPSILON = 1.0 / 1048576.0;

  function supertriangle(vertices) {
    var xmin = Number.POSITIVE_INFINITY,
        ymin = Number.POSITIVE_INFINITY,
        xmax = Number.NEGATIVE_INFINITY,
        ymax = Number.NEGATIVE_INFINITY,
        i, dx, dy, dmax, xmid, ymid;

    for(i = vertices.length; i--; ) {
      if(vertices[i][0] < xmin) xmin = vertices[i][0];
      if(vertices[i][0] > xmax) xmax = vertices[i][0];
      if(vertices[i][1] < ymin) ymin = vertices[i][1];
      if(vertices[i][1] > ymax) ymax = vertices[i][1];
    }

    dx = xmax - xmin;
    dy = ymax - ymin;
    dmax = Math.max(dx, dy);
    xmid = xmin + dx * 0.5;
    ymid = ymin + dy * 0.5;

    return [
      [xmid - 20 * dmax, ymid -      dmax],
      [xmid            , ymid + 20 * dmax],
      [xmid + 20 * dmax, ymid -      dmax]
    ];
  }

  function circumcircle(vertices, i, j, k) {
    var x1 = vertices[i][0],
        y1 = vertices[i][1],
        x2 = vertices[j][0],
        y2 = vertices[j][1],
        x3 = vertices[k][0],
        y3 = vertices[k][1],
        fabsy1y2 = Math.abs(y1 - y2),
        fabsy2y3 = Math.abs(y2 - y3),
        xc, yc, m1, m2, mx1, mx2, my1, my2, dx, dy;

    /* Check for coincident points */
    if(fabsy1y2 < EPSILON && fabsy2y3 < EPSILON)
      throw new Error("Eek! Coincident points!");

    if(fabsy1y2 < EPSILON) {
      m2  = -((x3 - x2) / (y3 - y2));
      mx2 = (x2 + x3) / 2.0;
      my2 = (y2 + y3) / 2.0;
      xc  = (x2 + x1) / 2.0;
      yc  = m2 * (xc - mx2) + my2;
    }

    else if(fabsy2y3 < EPSILON) {
      m1  = -((x2 - x1) / (y2 - y1));
      mx1 = (x1 + x2) / 2.0;
      my1 = (y1 + y2) / 2.0;
      xc  = (x3 + x2) / 2.0;
      yc  = m1 * (xc - mx1) + my1;
    }

    else {
      m1  = -((x2 - x1) / (y2 - y1));
      m2  = -((x3 - x2) / (y3 - y2));
      mx1 = (x1 + x2) / 2.0;
      mx2 = (x2 + x3) / 2.0;
      my1 = (y1 + y2) / 2.0;
      my2 = (y2 + y3) / 2.0;
      xc  = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
      yc  = (fabsy1y2 > fabsy2y3) ?
        m1 * (xc - mx1) + my1 :
        m2 * (xc - mx2) + my2;
    }

    dx = x2 - xc;
    dy = y2 - yc;
    return {i: i, j: j, k: k, x: xc, y: yc, r: dx * dx + dy * dy};
  }

  function dedup(edges) {
    var i, j, a, b, m, n;

    for(j = edges.length; j; ) {
      b = edges[--j];
      a = edges[--j];

      for(i = j; i; ) {
        n = edges[--i];
        m = edges[--i];

        if((a === m && b === n) || (a === n && b === m)) {
          edges.splice(j, 2);
          edges.splice(i, 2);
          break;
        }
      }
    }
  }

  Delaunay = {
    triangulate: function(vertices, key) {
      var n = vertices.length,
          i, j, indices, st, open, closed, edges, dx, dy, a, b, c;

      /* Bail if there aren't enough vertices to form any triangles. */
      if(n < 3)
        return [];

      /* Slice out the actual vertices from the passed objects. (Duplicate the
       * array even if we don't, though, since we need to make a supertriangle
       * later on!) */
      vertices = vertices.slice(0);

      if(key)
        for(i = n; i--; )
          vertices[i] = vertices[i][key];

      /* Make an array of indices into the vertex array, sorted by the
       * vertices' x-position. Force stable sorting by comparing indices if
       * the x-positions are equal. */
      indices = new Array(n);

      for(i = n; i--; )
        indices[i] = i;

      indices.sort(function(i, j) {
        var diff = vertices[j][0] - vertices[i][0];
        return diff !== 0 ? diff : i - j;
      });

      /* Next, find the vertices of the supertriangle (which contains all other
       * triangles), and append them onto the end of a (copy of) the vertex
       * array. */
      st = supertriangle(vertices);
      vertices.push(st[0], st[1], st[2]);
      
      /* Initialize the open list (containing the supertriangle and nothing
       * else) and the closed list (which is empty since we havn't processed
       * any triangles yet). */
      open   = [circumcircle(vertices, n + 0, n + 1, n + 2)];
      closed = [];
      edges  = [];

      /* Incrementally add each vertex to the mesh. */
      for(i = indices.length; i--; edges.length = 0) {
        c = indices[i];

        /* For each open triangle, check to see if the current point is
         * inside it's circumcircle. If it is, remove the triangle and add
         * it's edges to an edge list. */
        for(j = open.length; j--; ) {
          /* If this point is to the right of this triangle's circumcircle,
           * then this triangle should never get checked again. Remove it
           * from the open list, add it to the closed list, and skip. */
          dx = vertices[c][0] - open[j].x;
          if(dx > 0.0 && dx * dx > open[j].r) {
            closed.push(open[j]);
            open.splice(j, 1);
            continue;
          }

          /* If we're outside the circumcircle, skip this triangle. */
          dy = vertices[c][1] - open[j].y;
          if(dx * dx + dy * dy - open[j].r > EPSILON)
            continue;

          /* Remove the triangle and add it's edges to the edge list. */
          edges.push(
            open[j].i, open[j].j,
            open[j].j, open[j].k,
            open[j].k, open[j].i
          );
          open.splice(j, 1);
        }

        /* Remove any doubled edges. */
        dedup(edges);

        /* Add a new triangle for each edge. */
        for(j = edges.length; j; ) {
          b = edges[--j];
          a = edges[--j];
          open.push(circumcircle(vertices, a, b, c));
        }
      }

      /* Copy any remaining open triangles to the closed list, and then
       * remove any triangles that share a vertex with the supertriangle,
       * building a list of triplets that represent triangles. */
      for(i = open.length; i--; )
        closed.push(open[i]);
      open.length = 0;

      for(i = closed.length; i--; )
        if(closed[i].i < n && closed[i].j < n && closed[i].k < n)
          open.push(closed[i].i, closed[i].j, closed[i].k);

      /* Yay, we're done! */
      return open;
    },
    contains: function(tri, p) {
      /* Bounding box test first, for quick rejections. */
      if((p[0] < tri[0][0] && p[0] < tri[1][0] && p[0] < tri[2][0]) ||
         (p[0] > tri[0][0] && p[0] > tri[1][0] && p[0] > tri[2][0]) ||
         (p[1] < tri[0][1] && p[1] < tri[1][1] && p[1] < tri[2][1]) ||
         (p[1] > tri[0][1] && p[1] > tri[1][1] && p[1] > tri[2][1]))
        return null;

      var a = tri[1][0] - tri[0][0],
          b = tri[2][0] - tri[0][0],
          c = tri[1][1] - tri[0][1],
          d = tri[2][1] - tri[0][1],
          i = a * d - b * c;

      /* Degenerate tri. */
      if(i === 0.0)
        return null;

      var u = (d * (p[0] - tri[0][0]) - b * (p[1] - tri[0][1])) / i,
          v = (a * (p[1] - tri[0][1]) - c * (p[0] - tri[0][0])) / i;

      /* If we're outside the tri, fail. */
      if(u < 0.0 || v < 0.0 || (u + v) > 1.0)
        return null;

      return [u, v];
    }
  };

  if(typeof module !== "undefined")
    module.exports = Delaunay;
})();