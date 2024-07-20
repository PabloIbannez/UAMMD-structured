awk '{print $1,($1*20.62)/$3}' scaling2000.dat > batching2000.dat
awk '{print $1,($1*21.348709)/$3}' scaling5000.dat > batching5000.dat
awk '{print $1,($1*24.54)/$3}' scaling10000.dat > batching10000.dat

