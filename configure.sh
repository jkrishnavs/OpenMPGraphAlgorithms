#!/bin/bash
# if gcc version needed is < 6 TASKLOOP_DEFINED
rm Makefile.in
touch Makefile.in
if [ -z "$(gcc --version | head -n1 | cut -d" " -f4)"  ];
  then
      echo "Error: No gcc version found"
  elif [[ "$(gcc --version | head -n1 | cut -d" " -f4)" < "6.0.0"  ]];
  then
	echo "gcc version less than 6, Task loop might not work"
        echo "CCINST = gcc" > Makefile.in
	echo "TASKLOOP_DEFINED= No" >> Makefile.in
  else
      echo "GCC version is ","$(gcc --version | head -n1 | cut -d" " -f4)"
      echo "CCINST = gcc" > Makefile.in
      echo "TASKLOOP_DEFINED= yes" >> Makefile.in
fi

details=$(uname -a)
echo $details
# TODO: the below does not work
if [[ "$uname -a" == *"_64"* ]]; then
   echo CAPABILITY="GM_NODE64" >> Makefile.in
   else
   echo CAPABILITY="GM_NODE32" >> Makefile.in
fi


onlinecores=0
if [ $(ps -e -o psr= | sort | uniq | wc -l) -gt 32 ];
then
    echo "We might have problem in setting CPU FLAGS".
    
else
    onlinelist=$(lscpu | awk '/^On-line/{print $4}' | tail -1)
    corelist=$(echo $onlinelist | tr "," "\n")
    for coreset in $corelist
    do
	IFS='-' read -ra mybunch <<< "$coreset"
	if [ ${#mybunch[@]} == 1 ]; then
	    mycore=${mybunch[0]}
	    let "onlinecores= onlinecores + (1<<$mycore)"
	else
	    startcore=${mybunch[0]}
	    endcore=${mybunch[1]}
	    while [ $startcore -le $endcore ]
	    do
		let "onlinecores=onlinecores + (1 << $startcore)"
		let "startcore=startcore +1"
	    done
	fi
    done
fi
echo $onlinecores
echo "ONLINECORESFLAG = -D ONLINECORES=$onlinecores" >> Makefile.in

# Change the default chunk size here
echo "CHUNKSIZE = 1024" >> Makefile.in
echo "NUMTASKS = 1024" >> Makefile.in
echo "NUM_LOCKS = 1024" >> Makefile.in
