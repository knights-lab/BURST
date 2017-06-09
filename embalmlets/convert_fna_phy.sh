# depends on linelen in PATH
echo $(linelen $1 2) > $2
sed ':a;N;$!ba;s/\n>/#/g' $1 | sed ':a;N;$!ba;s/\n/ /g' | sed 's/#/\n/g' | sed 's/>//' >> $2

