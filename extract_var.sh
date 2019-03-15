#/bin/bash


exp_caso='LAMIP'

st=`echo $exp_caso | tr A-Z a-z`
ds=04
pathfile='/scratch/201803026n-2/cesm/'$exp_caso'/run/UNSET/atm/hist/*'
#pathscra00=/scratch/201803026n-2/cesm/var/$st/var00
#pathscra01=/scratch/201803026n-2/cesm/var/$st/var01
pathscra04=/scratch/201803026n-2/cesm/var/$st/var04

ifile=$pathfile

for f in $pathfile
do
    echo $f
    fstr=`echo  $f | tail -c 20`
    #cdo select,name=U,V $f $pathscra00/$exp_caso.var00.$fstr
    #cdo select,name=T $f $pathscra01/$exp_caso.var01.$fstr
    #cdo select,name=Q $f $pathscra02/$exp_caso.var02.$fstr
    cdo select,name=Z3 $f $pathscra04/$exp_caso.var04.$fstr
done
echo 'programa finalizado... '
