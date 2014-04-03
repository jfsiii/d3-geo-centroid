(function() {
    !function() {
        var d3 = {
            version: "3.4.4"
        };
        function d3_noop() {}
        var π = Math.PI, τ = 2 * π, halfπ = π / 2, ε = 1e-6, ε2 = ε * ε, d3_radians = π / 180, d3_degrees = 180 / π;
        function d3_sgn(x) {
            return x > 0 ? 1 : x < 0 ? -1 : 0;
        }
        function d3_cross2d(a, b, c) {
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        }
        function d3_acos(x) {
            return x > 1 ? 0 : x < -1 ? π : Math.acos(x);
        }
        function d3_asin(x) {
            return x > 1 ? halfπ : x < -1 ? -halfπ : Math.asin(x);
        }
        function d3_sinh(x) {
            return ((x = Math.exp(x)) - 1 / x) / 2;
        }
        function d3_cosh(x) {
            return ((x = Math.exp(x)) + 1 / x) / 2;
        }
        function d3_tanh(x) {
            return ((x = Math.exp(2 * x)) - 1) / (x + 1);
        }
        function d3_haversin(x) {
            return (x = Math.sin(x / 2)) * x;
        }
        d3.geo = {};
        d3.geo.stream = function(object, listener) {
            if (object && d3_geo_streamObjectType.hasOwnProperty(object.type)) {
                d3_geo_streamObjectType[object.type](object, listener);
            } else {
                d3_geo_streamGeometry(object, listener);
            }
        };
        function d3_geo_streamGeometry(geometry, listener) {
            if (geometry && d3_geo_streamGeometryType.hasOwnProperty(geometry.type)) {
                d3_geo_streamGeometryType[geometry.type](geometry, listener);
            }
        }
        var d3_geo_streamObjectType = {
            Feature: function(feature, listener) {
                d3_geo_streamGeometry(feature.geometry, listener);
            },
            FeatureCollection: function(object, listener) {
                var features = object.features, i = -1, n = features.length;
                while (++i < n) d3_geo_streamGeometry(features[i].geometry, listener);
            }
        };
        var d3_geo_streamGeometryType = {
            Sphere: function(object, listener) {
                listener.sphere();
            },
            Point: function(object, listener) {
                object = object.coordinates;
                listener.point(object[0], object[1], object[2]);
            },
            MultiPoint: function(object, listener) {
                var coordinates = object.coordinates, i = -1, n = coordinates.length;
                while (++i < n) object = coordinates[i], listener.point(object[0], object[1], object[2]);
            },
            LineString: function(object, listener) {
                d3_geo_streamLine(object.coordinates, listener, 0);
            },
            MultiLineString: function(object, listener) {
                var coordinates = object.coordinates, i = -1, n = coordinates.length;
                while (++i < n) d3_geo_streamLine(coordinates[i], listener, 0);
            },
            Polygon: function(object, listener) {
                d3_geo_streamPolygon(object.coordinates, listener);
            },
            MultiPolygon: function(object, listener) {
                var coordinates = object.coordinates, i = -1, n = coordinates.length;
                while (++i < n) d3_geo_streamPolygon(coordinates[i], listener);
            },
            GeometryCollection: function(object, listener) {
                var geometries = object.geometries, i = -1, n = geometries.length;
                while (++i < n) d3_geo_streamGeometry(geometries[i], listener);
            }
        };
        function d3_geo_streamLine(coordinates, listener, closed) {
            var i = -1, n = coordinates.length - closed, coordinate;
            listener.lineStart();
            while (++i < n) coordinate = coordinates[i], listener.point(coordinate[0], coordinate[1], coordinate[2]);
            listener.lineEnd();
        }
        function d3_geo_streamPolygon(coordinates, listener) {
            var i = -1, n = coordinates.length;
            listener.polygonStart();
            while (++i < n) d3_geo_streamLine(coordinates[i], listener, 1);
            listener.polygonEnd();
        }
        d3.geo.centroid = function(object) {
            d3_geo_centroidW0 = d3_geo_centroidW1 = d3_geo_centroidX0 = d3_geo_centroidY0 = d3_geo_centroidZ0 = d3_geo_centroidX1 = d3_geo_centroidY1 = d3_geo_centroidZ1 = d3_geo_centroidX2 = d3_geo_centroidY2 = d3_geo_centroidZ2 = 0;
            d3.geo.stream(object, d3_geo_centroid);
            var x = d3_geo_centroidX2, y = d3_geo_centroidY2, z = d3_geo_centroidZ2, m = x * x + y * y + z * z;
            if (m < ε2) {
                x = d3_geo_centroidX1, y = d3_geo_centroidY1, z = d3_geo_centroidZ1;
                if (d3_geo_centroidW1 < ε) x = d3_geo_centroidX0, y = d3_geo_centroidY0, z = d3_geo_centroidZ0;
                m = x * x + y * y + z * z;
                if (m < ε2) return [ NaN, NaN ];
            }
            return [ Math.atan2(y, x) * d3_degrees, d3_asin(z / Math.sqrt(m)) * d3_degrees ];
        };
        var d3_geo_centroidW0, d3_geo_centroidW1, d3_geo_centroidX0, d3_geo_centroidY0, d3_geo_centroidZ0, d3_geo_centroidX1, d3_geo_centroidY1, d3_geo_centroidZ1, d3_geo_centroidX2, d3_geo_centroidY2, d3_geo_centroidZ2;
        var d3_geo_centroid = {
            sphere: d3_noop,
            point: d3_geo_centroidPoint,
            lineStart: d3_geo_centroidLineStart,
            lineEnd: d3_geo_centroidLineEnd,
            polygonStart: function() {
                d3_geo_centroid.lineStart = d3_geo_centroidRingStart;
            },
            polygonEnd: function() {
                d3_geo_centroid.lineStart = d3_geo_centroidLineStart;
            }
        };
        function d3_geo_centroidPoint(λ, φ) {
            λ *= d3_radians;
            var cosφ = Math.cos(φ *= d3_radians);
            d3_geo_centroidPointXYZ(cosφ * Math.cos(λ), cosφ * Math.sin(λ), Math.sin(φ));
        }
        function d3_geo_centroidPointXYZ(x, y, z) {
            ++d3_geo_centroidW0;
            d3_geo_centroidX0 += (x - d3_geo_centroidX0) / d3_geo_centroidW0;
            d3_geo_centroidY0 += (y - d3_geo_centroidY0) / d3_geo_centroidW0;
            d3_geo_centroidZ0 += (z - d3_geo_centroidZ0) / d3_geo_centroidW0;
        }
        function d3_geo_centroidLineStart() {
            var x0, y0, z0;
            d3_geo_centroid.point = function(λ, φ) {
                λ *= d3_radians;
                var cosφ = Math.cos(φ *= d3_radians);
                x0 = cosφ * Math.cos(λ);
                y0 = cosφ * Math.sin(λ);
                z0 = Math.sin(φ);
                d3_geo_centroid.point = nextPoint;
                d3_geo_centroidPointXYZ(x0, y0, z0);
            };
            function nextPoint(λ, φ) {
                λ *= d3_radians;
                var cosφ = Math.cos(φ *= d3_radians), x = cosφ * Math.cos(λ), y = cosφ * Math.sin(λ), z = Math.sin(φ), w = Math.atan2(Math.sqrt((w = y0 * z - z0 * y) * w + (w = z0 * x - x0 * z) * w + (w = x0 * y - y0 * x) * w), x0 * x + y0 * y + z0 * z);
                d3_geo_centroidW1 += w;
                d3_geo_centroidX1 += w * (x0 + (x0 = x));
                d3_geo_centroidY1 += w * (y0 + (y0 = y));
                d3_geo_centroidZ1 += w * (z0 + (z0 = z));
                d3_geo_centroidPointXYZ(x0, y0, z0);
            }
        }
        function d3_geo_centroidLineEnd() {
            d3_geo_centroid.point = d3_geo_centroidPoint;
        }
        function d3_geo_centroidRingStart() {
            var λ00, φ00, x0, y0, z0;
            d3_geo_centroid.point = function(λ, φ) {
                λ00 = λ, φ00 = φ;
                d3_geo_centroid.point = nextPoint;
                λ *= d3_radians;
                var cosφ = Math.cos(φ *= d3_radians);
                x0 = cosφ * Math.cos(λ);
                y0 = cosφ * Math.sin(λ);
                z0 = Math.sin(φ);
                d3_geo_centroidPointXYZ(x0, y0, z0);
            };
            d3_geo_centroid.lineEnd = function() {
                nextPoint(λ00, φ00);
                d3_geo_centroid.lineEnd = d3_geo_centroidLineEnd;
                d3_geo_centroid.point = d3_geo_centroidPoint;
            };
            function nextPoint(λ, φ) {
                λ *= d3_radians;
                var cosφ = Math.cos(φ *= d3_radians), x = cosφ * Math.cos(λ), y = cosφ * Math.sin(λ), z = Math.sin(φ), cx = y0 * z - z0 * y, cy = z0 * x - x0 * z, cz = x0 * y - y0 * x, m = Math.sqrt(cx * cx + cy * cy + cz * cz), u = x0 * x + y0 * y + z0 * z, v = m && -d3_acos(u) / m, w = Math.atan2(m, u);
                d3_geo_centroidX2 += v * cx;
                d3_geo_centroidY2 += v * cy;
                d3_geo_centroidZ2 += v * cz;
                d3_geo_centroidW1 += w;
                d3_geo_centroidX1 += w * (x0 + (x0 = x));
                d3_geo_centroidY1 += w * (y0 + (y0 = y));
                d3_geo_centroidZ1 += w * (z0 + (z0 = z));
                d3_geo_centroidPointXYZ(x0, y0, z0);
            }
        }
        if (typeof define === "function" && define.amd) {
            define(d3);
        } else if (typeof module === "object" && module.exports) {
            module.exports = d3;
        } else {
            this.d3 = d3;
        }
    }();
})();