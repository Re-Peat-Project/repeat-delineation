import rasterio, numpy as np
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from _common import SRC_DTM as p, DATA
with rasterio.open(p) as s:
    for res in (10,30):
        f=res/2; h=int(s.height//f); w=int(np.ceil(s.width/f))
        a=s.read(1,out_shape=(h,w),resampling=Resampling.average)
        a[a==s.nodata]=np.nan
        t=from_origin(s.bounds.left,s.bounds.top,res,res)
        with rasterio.open(DATA / f"dtm_{res}m.tif","w",driver="GTiff",crs=s.crs,transform=t,width=w,height=h,count=1,dtype="float32",nodata=np.nan,compress="deflate") as d:
            d.write(a,1)
        print(res,a.shape,np.isfinite(a).mean(), np.nanmin(a), np.nanmax(a))
