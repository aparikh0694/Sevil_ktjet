#!/bin/sh

echo ' '
echo 'Backup :'
echo '========'
echo ' '

#time=`date -I`
time=`date "+%m-%d-%y%n"`
echo `date`

echo ' '
echo 'backup.sh'
echo ' '

#echo 'Documentation (doxygen & root) ...'
#doxygen doc/Doxyfile
#root -b -q 'make_html.C'
#ln -s htmldoc/ClassIndex.html htmldoc/index.html

tar czf Backup/$time.tar.gz examples Makefile configure *.h *.C *.sh 
#htmldoc

echo ''
echo 'Done and create file : Backup/'$time'.tar.gz' 
echo ' '
