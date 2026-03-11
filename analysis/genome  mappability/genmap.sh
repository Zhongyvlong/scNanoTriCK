genmap index -F hg38.fa -I genmap_index/hg38/

genmap map -K 150 -E 2 -I genmap_index/hg38/ \
	-O 150-kmer/hg38_mappability \
	-t -w -bg

genmap map -K 3000 -E 4 -I genmap_index/hg38/ \
	-O /3000-kmer/hg38_mappability \
	-t -w -bg
