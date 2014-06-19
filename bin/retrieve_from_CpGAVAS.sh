for ID in $(cut -f2 -d"=" $1)
do
	wget http://www.herbalgenomics.org/0506/cpgavas/docroot/tmp/dir_${ID}/${ID}.htm -r -l 1 -nH --cut-dir 5 -P $ID
done
