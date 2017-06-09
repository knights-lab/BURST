# process_ncbi.sh
# $1 is the input name, $2 is the output
# download genes from ftp or from nuccore (BioProject)
sed 's/; from TYPE material//g' $1 | sed 's/; from synonym TYPE material//g' | sed 's/; from verified material//g' | sed 's/; from reference material//g' | sed 's/:/-/g' | sed 's/ /_/g' | sed 's/,/_/g' | sed 's/#/_/g' | sed 's/\[/{/g' | sed 's/\]/}/g' | sed 's/;/_/g' | sed 's/\//_/g' | sed 's/</_/g' | sed -e '/^>/s/$/@/' -e 's/^>/#/'  | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort -u -t '	' -f -k 2,2  | sed -e 's/^/>/' -e 's/\t/\n/' | tail -n +2 > $2
