## Robust Rank Aggregation to measure group over representation

Example in run_RRA.sh:
```
#!/bin/bash

location_of_tool="RRA_tool.py"
location_input="test_input.tsv"
columnsInGroup="test_group1.txt"
columnsOutOfGroup="test_group2.txt"
enrichedInHighOrLowValues="high"
outputPrefix="test_output"

python2.7 ${location_of_tool} \
--location_input ${location_input} \
--columnsInGroup ${columnsInGroup} \
--columnsOutOfGroup ${columnsOutOfGroup} \
--enrichedInHighOrLowValues ${enrichedInHighOrLowValues} \
--outputPrefix ${outputPrefix}
```
