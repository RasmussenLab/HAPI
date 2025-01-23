# Author: Kirstine Ravn
for RSID in rs333 rs143241023 rs150628438 rs369842709 rs556322139 rs61231801 rs66552573 rs67580019;
do
    echo "############# Running $RSID #############"
    script_name="PyScriptDelres_no_duplicated_RISE/Py${RSID}/Py${RSID}.py" 
    python3 $script_name --deletion_rsid $RSID --path PyScriptDelres_no_duplicated_RISE
    echo "############# Done running $RSID #############"
done




for RSID in rs333 rs143241023 rs150628438 rs369842709 rs556322139 rs61231801 rs66552573 rs67580019;
do
    echo "Running $RSID"
    python3 process_hapi_results.py --deletion_rsid $RSID --path PyScriptDelres_no_duplicated_RISE
done
