The files in this directory are derived from the [Genome Taxonomy Database](https://gtdb.ecogenomic.org/) version r220, where they are shared under the Creative Commons Attribution-ShareAlike 4.0 International License. Additional information can be found at https://gtdb.ecogenomic.org/licenses.

## Code used to derive the example data

```python
import skbio
from random import choices
```

```python
l = list(skbio.io.read(open('./ssu_all_r220.fna'), format='fasta'))
l = [e for e in l if not skbio.DNA(e).has_degenerates()]
```

```python
query_n = 5

query = open(f'query.fasta', 'w')
for e in choices(l, k=query_n):
    e[:250].write(query, format='fasta')
query.close()
```

```python
reference_n = 20

reference = open('reference.fasta', 'w')
reference_md = open('reference-metadata.tsv', 'w')
reference_md.write('id\tdescription\n')
for e in choices(l, k=reference_n):
    e.write(reference, format='fasta')
    reference_md.write(f"{e.metadata['id']}\t"
                        f"{e.metadata['description']}\n")
reference_md.close()
reference.close()
```

