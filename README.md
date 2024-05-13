# CNVizard
A streamlit App for the analysis of copy number variants

# Setup (mandatory)
1. Create a virtual enviroment and activate it
2. Navigate into the CNVizard folder 
3. Install the requirements : pip install requirements.txt
4. Start the app : streamlit run run_visualizer.py

# Setup (optional) 
1. Create an .env file in the CNVizard folder (optional). This can be used to set up IGV outlinks : APPSETTING_IGV_OUTLINK = "/path/to/samplename/samplename.cram" (samplename needs to be included in the path in any way or form)
2. Change the AnnotSV table format : Navigate to /CNVizard/Ressources/Annotsv_format.txt -> Change the standard column Names
3. Add/Remove a panel-list : Navigate to /CNVizard/Ressources/Candi_lists/ -> Add/Remove new panel.txt file
4. Modiy OMIM-List : Navigate to /CNVizard/Ressources/OMIM/omim.txt -> Modify .txt file
5. Change reference List : Navigate tot /CNVizard/Ressources/References/ -> replace bintest_ref.parquet and/or total_ref.parquet
