# run with
# python *.py allpaths.txt capitalist-lca.txt taxonomymap.txt bt2out.sam outfile.txt
import sys
from collections import defaultdict

# build sets of tied best hits for each query
hitsets = defaultdict(set)
hitid = {}

print('loading query-ref map')
with open(sys.argv[1],'r') as f:
    for line in f:
        words = line.strip().split('\t')
        hitsets[words[0]].add(words[1])
        hitid[words[0]] = words[2]


print('loading query-tax LCA map')
hit2tax = {}
with open(sys.argv[2],'r') as f:
    for line in f:
        words = line.strip().split('\t')
        hit2tax[words[0]] = words[12]

print('loading ref2tax map')
ref2tax = {}
with open(sys.argv[3],'r') as f:
    for line in f:
        words = line.strip().split('\t')
        ref2tax[words[0]] = words[1]

print('parsing bt2 output')
hittype = {}
hittaxtype = {}

# track whether each hit was found
for f in hitsets:
    hittype[f] = 'Miss'

for f in hitsets:
    hittaxtype[f] = 'Miss'

with open(sys.argv[4],'r') as f:
    for line in f:
        words = line.strip().split('\t')
        query = words[0]
        ref = words[2]
        if ref == '*':
            continue
        truetax = hit2tax[query]
        obstax = ref2tax[ref]
        if ref in hitsets[query]:
            hittype[query] = 'Correct'
        else:
            hittype[query] = 'Incorrect'

        if obstax == truetax:
            hittaxtype[query] = 'Correct'
        elif truetax in obstax:
            hittaxtype[query] = 'Over-specific'
        else:
            hittaxtype[query] = 'Incorrect'
            
print('Writing output file...')
with open(sys.argv[5],'w') as f:
    f.write('Query\tOptimal ID\tMatch\tTaxonomy\n')
    for hit in hittype:
        f.write(hit + '\t' + hitid[hit] + '\t' + hittype[hit] + '\t' + hittaxtype[hit] + '\n')

