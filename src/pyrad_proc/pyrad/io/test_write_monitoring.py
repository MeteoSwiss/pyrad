#%%
import pyrad
import datetime
pyrad.io.write_monitoring_ts(datetime.datetime.now(), 2000, [1,2,3], [10,20,30], "Test", "/users/wolfensb/test.csv")

