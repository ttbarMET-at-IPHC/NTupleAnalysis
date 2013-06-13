#!/bin/bash

user=`whoami`

echo "#########################################"
echo "  Cleaning PROOF instances"
echo "#########################################"
echo "Seach process to be killed "
ps -f -u $user | grep proofserv.exe | awk '{print "kill -9 " $2}' > kill.sh
chmod +x kill.sh
echo "Kill them"
./kill.sh
echo "Done"
echo "#########################################"

