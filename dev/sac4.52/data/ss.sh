#!/bin/sh

DIR_WORK=`pwd`
# nujno vvesti kolichestvo picture
n=$1
i=1  
g=1 # step
  
  cd ../../../
  echo $pwd
  echo ------------------;
    while [ $i -le $n  ]
    do
        #echo "$n * $i = `expr $i \* $n`"   
	t=`printf %05d $i`
	
	convert DATA/opoz16/$t.png /data/ap1vf/gif/opoz16/$t.gif
	i=`expr $i + $g`
      echo "$t"   	
    done 

    cd DATA/gif/$b
    gifsicle --delay 50 --loop --optimize *.gif > test_$b.gif
    cd $DIR_WORK
