use crate::inc::inc::BigFloatInc;

// arcsine polynomial coefficients
pub(crate) const ASIN_VALUES: [BigFloatInc; 20] = [
    BigFloatInc {
        sign: -1,
        e: -44,
        n: 44,
        m: [3932, 2141, 8173, 5080, 603, 3474, 2073, 7579, 197, 5113, 1178],
    },
    BigFloatInc {
        sign: -1,
        e: -45,
        n: 44,
        m: [1348, 7318, 8389, 8931, 6357, 316, 2165, 9553, 2944, 6504, 2651],
    },
    BigFloatInc {
        sign: -1,
        e: -46,
        n: 44,
        m: [5439, 7732, 3302, 7773, 5112, 8799, 4776, 4146, 5431, 8167, 7891],
    },
    BigFloatInc {
        sign: -1,
        e: -46,
        n: 44,
        m: [7684, 1936, 7582, 4311, 2781, 7994, 6208, 7452, 6778, 4098, 2685],
    },
    BigFloatInc {
        sign: -1,
        e: -47,
        n: 44,
        m: [5565, 4403, 4734, 4511, 7513, 2160, 8314, 3802, 7685, 1908, 9887],
    },
    BigFloatInc {
        sign: -1,
        e: -47,
        n: 44,
        m: [6357, 361, 6163, 4858, 8298, 4491, 6365, 5737, 3621, 4554, 3834],
    },
    BigFloatInc {
        sign: -1,
        e: -47,
        n: 44,
        m: [9915, 9835, 8717, 97, 9982, 5735, 7156, 6808, 3028, 9118, 1542],
    },
    BigFloatInc {
        sign: -1,
        e: -48,
        n: 44,
        m: [9501, 1306, 7418, 6507, 2572, 6224, 1658, 5955, 982, 5287, 6381],
    },
    BigFloatInc {
        sign: -1,
        e: -48,
        n: 44,
        m: [9424, 3446, 1131, 5615, 3718, 6109, 861, 4826, 7710, 2891, 2696],
    },
    BigFloatInc {
        sign: -1,
        e: -48,
        n: 44,
        m: [1264, 1136, 3998, 8, 8610, 8589, 8298, 1478, 7254, 7623, 1158],
    },
    BigFloatInc {
        sign: -1,
        e: -49,
        n: 44,
        m: [5786, 9061, 8647, 3514, 1729, 2926, 4306, 2017, 9299, 5474, 5049],
    },
    BigFloatInc {
        sign: -1,
        e: -49,
        n: 44,
        m: [9792, 2852, 1362, 6891, 9545, 239, 6840, 647, 5316, 88, 2226],
    },
    BigFloatInc {
        sign: -1,
        e: -50,
        n: 44,
        m: [1968, 1361, 2973, 8696, 2934, 726, 3312, 5390, 4462, 2274, 9909],
    },
    BigFloatInc {
        sign: -1,
        e: -50,
        n: 44,
        m: [1924, 6750, 4930, 9950, 1846, 6582, 193, 4230, 1621, 1692, 4448],
    },
    BigFloatInc {
        sign: -1,
        e: -50,
        n: 44,
        m: [3074, 6869, 347, 9155, 9017, 8529, 248, 90, 260, 2421, 2011],
    },
    BigFloatInc {
        sign: -1,
        e: -51,
        n: 44,
        m: [9376, 6322, 2178, 1477, 1838, 1464, 9939, 5295, 8389, 5324, 9151],
    },
    BigFloatInc {
        sign: -1,
        e: -51,
        n: 44,
        m: [8484, 3657, 5757, 6314, 1748, 716, 2762, 9486, 8863, 4028, 4187],
    },
    BigFloatInc {
        sign: -1,
        e: -51,
        n: 44,
        m: [8319, 8757, 3585, 4600, 4182, 8925, 4565, 4448, 1568, 5137, 1925],
    },
    BigFloatInc {
        sign: -1,
        e: -52,
        n: 44,
        m: [7972, 2629, 2365, 3272, 8447, 9098, 2073, 5229, 1515, 4827, 8893],
    },
    BigFloatInc {
        sign: -1,
        e: -52,
        n: 44,
        m: [6102, 4039, 8633, 8733, 9459, 1758, 5928, 7531, 6638, 814, 4124],
    },
];
