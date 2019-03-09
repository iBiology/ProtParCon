for t in `ls test_*.py`
do
    echo "Testing $t"
    python $t
done