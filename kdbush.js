(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
typeof define === 'function' && define.amd ? define(factory) :
(global.KDBush = factory());
}(this, (function () { 'use strict';

function sortKD(ids, coords, nodeSize, left, right, axis, axisCount) {

    if (right - left <= nodeSize) {
        return;
    }

    var m = (left + right) >> 1; // middle index

    // sort ids and coords around the middle index so that the halves lie
    // either left/right or top/bottom correspondingly (taking turns)
    select(ids, coords, m, left, right, axis, axisCount);

    // recursively kd-sort first half and second half on the opposite axis
    sortKD(ids, coords, nodeSize, left, m - 1, (1 + axis) % axisCount, axisCount);
    sortKD(ids, coords, nodeSize, m + 1, right, (1 + axis) % axisCount, axisCount);
}

// custom Floyd-Rivest selection algorithm: sort ids and coords so that
// [left..k-1] items are smaller than k-th item (on either x or y axis)
function select(ids, coords, k, left, right, axis, axisCount) {

    while (right > left) {
        if (right - left > 600) {
            var n = right - left + 1;
            var m = k - left + 1;
            var z = Math.log(n);
            var s = 0.5 * Math.exp(2 * z / 3);
            var sd = 0.5 * Math.sqrt(z * s * (n - s) / n) * (m - n / 2 < 0 ? -1 : 1);

            var newLeft = Math.max(left, Math.floor(k - m * s / n + sd));
            var newRight = Math.min(right, Math.floor(k + (n - m) * s / n + sd));

            select(ids, coords, k, newLeft, newRight, axis, axisCount);
        }

        var t = coords[axisCount * k + axis];
        var i = left;
        var j = right;

        swapItem(ids, coords, left, k, axisCount);

        if (coords[axisCount * right + axis] > t) {
            swapItem(ids, coords, left, right, axisCount);
        }

        while (i < j) {
            swapItem(ids, coords, i, j, axisCount);
            i++;
            j--;
            while (coords[axisCount * i + axis] < t) { i++; }
            while (coords[axisCount * j + axis] > t) { j--; }
        }

        if (coords[axisCount * left + axis] === t) {
            swapItem(ids, coords, left, j, axisCount);
        } else {
            j++;
            swapItem(ids, coords, j, right, axisCount);
        }

        if (j <= k) { left = j + 1; }
        if (k <= j) { right = j - 1; }
    }
}

function swapItem(ids, coords, i, j, axisCount) {
    swap(ids, i, j);
    swap(coords, axisCount * i + 0, axisCount * j + 0);
    swap(coords, axisCount * i + 1, axisCount * j + 1);
    swap(coords, axisCount * i + 2, axisCount * j + 2);
}

function swap(arr, i, j) {
    var tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

function range(ids, coords, minX, minY, minZ, maxX, maxY, maxZ, nodeSize, axisCount) {
    var stack = [0, ids.length - 1, 0];
    var result = [];

    // recursively search for items in range in the kd-sorted arrays
    while (stack.length) {
        var axis = stack.pop();
        var right = stack.pop();
        var left = stack.pop();

        // if we reached "tree node", search linearly
        if (right - left <= nodeSize) {
            for (var i = left; i <= right; i++) {
                var x = coords[axisCount * i + 0];
                var y = coords[axisCount * i + 1];
                var z = coords[axisCount * i + 2];
                if (
                    x >= minX && x <= maxX &&
                    y >= minY && y <= maxY &&
                    z >= minZ && z <= maxZ
                ) { result.push(ids[i]); }
            }
            continue;
        }

        // otherwise find the middle index
        var m = (left + right) >> 1;

        // include the middle item if it's in range
        var x$1 = coords[axisCount * m + 0];
        var y$1 = coords[axisCount * m + 1];
        var z$1 = coords[axisCount * m + 2];
        if (
            x$1 >= minX && x$1 <= maxX &&
            y$1 >= minY && y$1 <= maxY &&
            z$1 >= minZ && z$1 <= maxZ
    ) { result.push(ids[m]); }

        // queue search in halves that intersect the query
        var next_axis = (1 + axis) % axisCount;
        var min_conditional = (void 0);
        var max_conditional = (void 0);

        switch (axis) {

            case 0:
                min_conditional = minX <= x$1;
                max_conditional = maxX >= x$1;
                break;

            case 1:
                min_conditional = minY <= y$1;
                max_conditional = maxY >= y$1;
                break;

            case 2:
                min_conditional = minZ <= z$1;
                max_conditional = maxZ >= z$1;
                break;

        }

        if (min_conditional) {
            stack.push(left);
            stack.push(m - 1);
            stack.push(next_axis);
        }

        if (max_conditional) {
            stack.push(m + 1);
            stack.push(right);
            stack.push(next_axis);
        }
    }

    return result;
}

function within(ids, coords, qx, qy, qz, r, nodeSize, axisCount) {
    var stack = [0, ids.length - 1, 0];
    var result = [];
    var r2 = r * r;

    // recursively search for items within radius in the kd-sorted arrays
    while (stack.length) {
        var axis = stack.pop();
        var right = stack.pop();
        var left = stack.pop();

        // if we reached "tree node", search linearly
        if (right - left <= nodeSize) {
            for (var i = left; i <= right; i++) {
                if (sqDist(
                    coords[axisCount * i + 0],
                    coords[axisCount * i + 1],
                    coords[axisCount * i + 2],
                    qx,
                    qy,
                    qz
                ) <= r2) {
                    result.push(ids[i]);
                }
            }
            continue;
        }

        // otherwise find the middle index
        var m = (left + right) >> 1;

        // include the middle item if it's in range
        var x = coords[axisCount * m + 0];
        var y = coords[axisCount * m + 1];
        var z = coords[axisCount * m + 2];
        if (sqDist(
            x,
            y,
            z,
            qx,
            qy,
            qz
        ) <= r2) {
            result.push(ids[m]);
        }

        // queue search in halves that intersect the query
        var next_axis = (1 + axis) % axisCount;
        var min_conditional = (void 0);
        var max_conditional = (void 0);

        switch (axis) {

            case 0:
                min_conditional = qx - r <= x;
                max_conditional = qx + r >= x;
                break;

            case 1:
                min_conditional = qy - r <= y;
                max_conditional = qy + r >= y;
                break;

            case 2:
                min_conditional = qz - r <= z;
                max_conditional = qz + r >= z;
                break;

        }

        if (min_conditional) {
            stack.push(left);
            stack.push(m - 1);
            stack.push(next_axis);
        }

        if (max_conditional) {
            stack.push(m + 1);
            stack.push(right);
            stack.push(next_axis);
        }
    }

    return result;
}

function sqDist(ax, ay, az, bx, by, bz) {
    var dx = ax - bx;
    var dy = ay - by;
    var dz = az - bz;
    return dx * dx + dy * dy + dz * dz;
}

var KDBush = function KDBush(ref) {
    var points = ref.points;
    var getX = ref.getX;
    var getY = ref.getY;
    var getZ = ref.getZ;
    var nodeSize = ref.nodeSize;
    var ArrayType = ref.ArrayType;
    var axisCount = ref.axisCount;

    this.nodeSize = nodeSize;
    this.points = points;
    this.axisCount = axisCount;

    var IndexArrayType = points.length < 65536 ? Uint16Array : Uint32Array;

    // store indices to the input array and coordinates in separate typed arrays
    var ids = this.ids = new IndexArrayType(points.length);
    var coords = this.coords = new ArrayType(points.length * axisCount);

    for (var i = 0; i < points.length; i++) {
        ids[i] = i;
        coords[axisCount * i + 0] = getX(points[i]);
        coords[axisCount * i + 1] = getY(points[i]);
        coords[axisCount * i + 2] = getZ(points[i]);
    }

    // kd-sort both arrays for efficient search (see comments in sort.js)
    sortKD(ids, coords, nodeSize, 0, ids.length - 1, 0, axisCount);

};

KDBush.prototype.range = function range$1 (minX, minY, minZ, maxX, maxY, maxZ) {
    return range(this.ids, this.coords, minX, minY, minZ, maxX, maxY, maxZ, this.nodeSize, this.axisCount);
};

KDBush.prototype.within = function within$1 (x, y, z, r) {
    return within(this.ids, this.coords, x, y, z, r, this.nodeSize, this.axisCount);
};

return KDBush;

})));
