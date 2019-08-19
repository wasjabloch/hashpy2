# Script to plot all possible fault mechanisms returned by HASH.

hypf=$1     #  NonLinLoc .hyp file
allfpsf=$2  #  All faulplane solutions
angf=$3     #  Take-off angles
outfile=$4     #  Postscriptfile 

J="A0/0/10c"
F="+f12p,Helvetica+jLT"
bbcolor="gray"

name=`basename $angf | cut -d. -f1`
T0=`grep GEOGRAPHIC $hypf | awk '{printf("%04d-%02d-%02d %02d:%02d:%02d",$3, $4, $5, $6, $7, $6, $7, $8)}'`
lat=`grep GEOGRAPHIC $hypf | awk '{printf("%.3f\260N", $10)}'`
lon=`grep GEOGRAPHIC $hypf | awk '{printf("%.3f\260E", $12)}'`
dep=`grep GEOGRAPHIC $hypf | awk '{printf("%.0f km", $14)}'`
sdr=`grep FOCALMECH $hypf | awk '{printf("%d %d %d", $6, $7, $8)}'`
rms=`grep FOCALMECH $hypf | awk '{printf("%.0f", $10)}'`

gmt psxy -Rg -J$J -K -T -Y10c > $outfile

#  plot prefered solution
echo $sdr | awk '{print 0, 0, 0, $1, $2, $3, 5, 0, 0}' |
  gmt psmeca -J$J -Rg -Sa10c -G$bbcolor -W0 -K -O  >> $outfile


#  plot all solutions
cat $allfpsf | awk '{print 0, 0, 0, $1, $2, $3, 5, 0, 0}' |
  gmt psmeca -J$J -Rg -Sa10c -Wthinner,100 -t50 -T -K -O >> $outfile


#  plot take off angles
awk '{
bearing=$2; plunge=$3-90;
if(plunge<0){plunge=-plunge; bearing+=180}
if(bearing<0){bearing=bearing+360;}
if(bearing>360){bearing=bearing-360;}
printf("%i/%i\n",plunge, bearing);}' $angf |
python stereonet.py --lines |
  psxy -J$J -Rd -Sp4p -Gsalmon -K -O >> $outfile


#  get coordinates
grep best_angle $angf | awk '{
bearing=$2; plunge=$3-90;
if(plunge<0){plunge=-plunge; bearing+=180}
if(bearing<0){bearing=bearing+360;}
if(bearing>360){bearing=bearing-360;}
printf("%i/%i\n",plunge, bearing);}'|
python stereonet.py --lines > /tmp/coordinates


#  get station names
grep best_angle $angf | awk '{print $1}' > /tmp/names


#  get polarities
grep best_angle $angf | awk '{
if($4==-1) s="-";
if($4==1) s="+";
print 0.3, s}' > /tmp/polarities


#  get p/s ratios
grep best_angle $angf | awk '{
if($5==0){print 0, "c"}
else {print 2*sqrt(1/$5), "c";}}' > /tmp/rings


#  Write Rings
paste /tmp/coordinates /tmp/rings |
  psxy -J$J -Rd -S -Wthick -K -O >> $outfile
#  Write Polarities
paste /tmp/coordinates /tmp/polarities |
  psxy -J$J -Rd -S -Wthick -K -O >> $outfile
#  Write Station names
paste /tmp/coordinates /tmp/names |
  pstext -J$J -Rd -D0.2/-0.2 -F+f8p,Helvetica,darkblue+jTL -N -K -O >> $outfile


#  get axis
echo $sdr | python3 strdiprake2ptnaxes.py -l | awk '{print $1}' |  
  python stereonet.py --lines | awk '{print $1, $2, "P"}' |
   pstext -J$J -Rd -F+f12p,Helvetica-Bold -N -K -O >> $outfile

echo $sdr | python3 strdiprake2ptnaxes.py -l | awk '{print $2}' |  
  python stereonet.py --lines | awk '{print $1, $2, "T"}' |
  pstext -J$J -Rd -F+f12p,Helvetica-Bold -N -K -O >> $outfile


#  Plot Infos
echo "0 0 $T0" | gmt pstext -J$J -Rg -F$F -N -K -O -Y2.5c -X6c >> $outfile
echo "0 0 ID: $name" | gmt pstext -J$J -Rg -F$F -N -K -O -Y-14p  >> $outfile
echo "0 0 Lat.: $lat" | gmt pstext -J$J -Rg -F$F -N -K -O -Y-14p  >> $outfile
echo "0 0 Lon.: $lon" | gmt pstext -J$J -Rg -F$F -N -K -O -Y-14p  >> $outfile
echo "0 0 Dep.: $dep" | gmt pstext -J$J -Rg -F$F -N -K -O -Y-14p  >> $outfile
echo "0 0 Str, Dip, Rake: $sdr" | gmt pstext -J$J -Rg -F$F -N -K -O -Y-14p  >> $outfile
echo "0 0 RMS: $rms" | gmt pstext -J$J -Rg -F$F -N -K -O -Y-14p  >> $outfile

gmt psxy -Rg -J$J -T -O >> $outfile

