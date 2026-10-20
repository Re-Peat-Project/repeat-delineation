import rasterio, numpy as np
from _common import SRC_DTM as p
with rasterio.open(p) as s:
    print(s.crs, s.res, s.shape, s.dtypes, s.nodata, s.bounds, s.overviews(1), s.profile.get('compress'), s.block_shapes)
    a=s.read(1,out_shape=(s.height//50,s.width//50))
    print(a.shape, np.isnan(a).mean() if a.dtype.kind=='f' else None)
    v=a[np.isfinite(a)&(a!=s.nodata)] if s.nodata is not None else a[np.isfinite(a)]
    print(len(v)/a.size, np.percentile(v,[0,1,5,25,50,75,95,99,100]))
    print((a==s.nodata).mean() if s.nodata is not None else 'no nodata')
