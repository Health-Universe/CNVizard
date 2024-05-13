"""
File which contains export functions of the CNVizard
@author: Jeremias Krause / Matthias Begemann / Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""
import pandas as pd 
import xlsxwriter
from io import BytesIO
class cnv_exporter:
    """
    Class which is used to export the filtered dataframe as excel tables
    """

    def __init__(self):
        """
        Constructor of the class cnv_exporter
        """

    def save_tables_as_excel(self,total_df:pd.DataFrame,bintest_df:pd.DataFrame,hom_del_df:pd.DataFrame,total_candi_df:pd.DataFrame,
                            bintest_candi_df:pd.DataFrame,consecutive_del_df:pd.DataFrame, consecutive_dup_df:pd.DataFrame):
        """
        Function which writes the preset DataFrames to an excel file 
        each preset represents an excel sheet
        Input Arguments: 
            total_df (pandas DataFrame) -> .cnr DataFrame 
            bintest_df (pandas DataFrame) -> bintest DataFrame
            hom_del_df (pandas DataFrame) -> .cnr DataFrame filtered for homozygous deletions
            total_candi_df (pandas DataFrame) -> .cnr DataFrame filtered for candigenes 
            bintest_candi_df (pandas DataFrame) -> bintest DataFrame filtered for candigenes
            consecutive_del_df (pandas DataFrame) -> .cnr DataFrame filtered for consecutive deletions
            consecutive_dup_df (pandas DataFrame) -> .cnr DataFrame filtered for consecutive duplications
        Output Arguments:
            processed_data (Output extracted from Output-Stream) -> Output Data which will be passed to the streamlit download Button object
        """
        #Open Output-Stream 
        output=BytesIO()
        #Write Data from Output-Stream to Excel file 
        writer = pd.ExcelWriter(output,engine="xlsxwriter")
        #Declare a list which contains the preset names which will be used to name the individual excel sheets
        list_of_saved_results =["total","bintest","hom_del","total_candi","bintest_candi","consecutive_del","consecutive_dup"]
        #Declare a list which contains the preset dataframes
        list_of_selections =[total_df,bintest_df,hom_del_df,total_candi_df,bintest_candi_df,consecutive_del_df,consecutive_dup_df]
        #Iterate over the sheetnames and DataFrames, format them and pass them to the writer
        for name,df in zip(list_of_saved_results,list_of_selections):
            #Open Excel writer
            df.to_excel(writer, sheet_name=name, header=False, index=False, startrow=1)
            #get number of columns of DataFrame
            num_columns = len(df.columns.tolist())
            #Create writer book and select current Sheet
            workbook = writer.book 
            worksheet = writer.sheets[name]
            #Add Formate Options to workbook
            header_format = workbook.add_format(
                {
                "bold": True,
                "text_wrap": True,
                "valign": "top",
                "fg_color": "#D7E4BC",
                "border": 1,
                }
            )
            #Add conditional color coding, similar to the pandas styler
            format_yellow = workbook.add_format({"bg_color": "#FFFF00"})
            format_red = workbook.add_format({"bg_color": "#ff0000"})
            worksheet.set_column(1, num_columns, 11)
            worksheet.conditional_format('G1:G1048576', {'type': '3_color_scale', 'min_type': 'num', 'min_value': 0,'mid_type': 'num','mid_value': 0.5,'max_type': 'num','max_value': 1})
            worksheet.conditional_format('I1:I1048576', {'type': 'cell', 'criteria': 'less than or equal to', 'value':  -0.65, 'format':   format_yellow})
            column_length = str(len(df))
            area_to_color = 'F1:' + "F" + column_length 
            worksheet.conditional_format(area_to_color, {'type': 'cell', 'criteria': 'equal to', 'value':  0, 'format':   format_red})
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)
            #Add Background colors
            bg_format1 = workbook.add_format({'bg_color': '#EEEEEE'}) # grey cell background color
            bg_format2 = workbook.add_format({'bg_color': '#FFFFFF'}) # white cell background color 
            for rn in range(len(df)):
                worksheet.set_row(rn, cell_format=(bg_format1 if rn%2==0 else bg_format2))
        #Close writer
        writer.close()
        #Extract output
        processed_data=output.getvalue()
        return processed_data

    def save_filtered_table_as_excel(self,filtered_df:pd.DataFrame, filtered_name:str):
        """
        Function which writes the preset DataFrames to an excel file 
        each preset represents an excel sheet
        Input Arguments: 
            filtered_df (pandas DataFrame) -> cnv dataframe which was filtered with the filter criteria provided by the user 
            filtered_name (String) -> name which will be used to name the excel sheet
        Output Arguments:
            processed_data (Output extracted from Output-Stream) -> Output Data which will be passed to the streamlit download Button object
        """
        #Open Output-Stream 
        output=BytesIO()
        #Write Data from Output-Stream to Excel file 
        writer = pd.ExcelWriter(output,engine="xlsxwriter")
        #Declare a list which contains the preset names which will be used to name the individual excel sheets
        list_of_saved_results =[filtered_name]
        #Declare a list which contains the preset dataframes
        list_of_selections =[filtered_df]
        #Iterate over the sheetnames and DataFrames, format them and pass them to the writer
        for name,df in zip(list_of_saved_results,list_of_selections):
            #Open Excel writer
            df.to_excel(writer, sheet_name=name, header=False, index=False, startrow=1)
            #get number of columns of DataFrame
            num_columns = len(df.columns.tolist())
            #Create writer book and select current Sheet
            workbook = writer.book 
            worksheet = writer.sheets[name]
            #Add Formate Options to workbook
            header_format = workbook.add_format(
                {
                "bold": True,
                "text_wrap": True,
                "valign": "top",
                "fg_color": "#D7E4BC",
                "border": 1,
                }
            )
            #Add conditional color coding, similar to the pandas styler
            format_yellow = workbook.add_format({"bg_color": "#FFFF00"})
            format_red = workbook.add_format({"bg_color": "#ff0000"})
            worksheet.set_column(1, num_columns, 11)
            worksheet.conditional_format('G1:G1048576', {'type': '3_color_scale', 'min_type': 'num', 'min_value': 0,'mid_type': 'num','mid_value': 0.5,'max_type': 'num','max_value': 1})
            worksheet.conditional_format('I1:I1048576', {'type': 'cell', 'criteria': 'less than or equal to', 'value':  -0.65, 'format':   format_yellow})
            column_length = str(len(df))
            area_to_color = 'F1:' + "F" + column_length 
            worksheet.conditional_format(area_to_color, {'type': 'cell', 'criteria': 'equal to', 'value':  0, 'format':   format_red})
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)
            #Add Background colors
            bg_format1 = workbook.add_format({'bg_color': '#EEEEEE'}) # grey cell background color
            bg_format2 = workbook.add_format({'bg_color': '#FFFFFF'}) # white cell background color 
            for rn in range(len(df)):
                worksheet.set_row(rn, cell_format=(bg_format1 if rn%2==0 else bg_format2))
        #Close writer
        writer.close()
        #Extract output
        processed_data=output.getvalue()
        return processed_data

    def save_tables_as_excel_tsv(self,filtered_tsv:pd.DataFrame):
        """
        slightly altered save_tables_as_excel function, used to export .tsv file 
        Function which writes the preset DataFrames to an excel file 
        each preset represents an excel sheet
        Input Arguments:
            filtered_tsv (pandas DataFrame) -> .tsv file annotated by AnnotSV, subsequently filtered and formatted
        Output Arguments:
            processed_data (Output extracted from Output-Stream) -> Output Data which will be passed to the streamlit download Button object
        """
        #Open Output-Stream
        output=BytesIO()
        #Write Data from Output-Stream to Excel file
        writer = pd.ExcelWriter(output,engine="xlsxwriter")
        #Declare a list which contains the preset names which will be used to name the individual excel sheets
        list_of_saved_results =["filered"]
        #Declare a list which contains the preset dataframes
        list_of_selections =[filtered_tsv]
        #Iterate over the sheetnames and DataFrames, format them and pass them to the writer
        for name,df in zip(list_of_saved_results,list_of_selections):
            #Open Excel writer
            df.to_excel(writer, sheet_name=name, header=False, index=False, startrow=1)
            #get number of columns of DataFrame
            num_columns = len(df.columns.tolist())
            #Create writer book and select current Sheet
            workbook = writer.book 
            worksheet = writer.sheets[name]
            #Add Formate Options to workbook
            header_format = workbook.add_format(
                {
                "bold": True,
                "text_wrap": True,
                "valign": "top",
                "fg_color": "#D7E4BC",
                "border": 1,
                }
            )
            #Add conditional color coding, similar to the pandas styler
            format_yellow = workbook.add_format({"bg_color": "#FFFF00"})
            format_red = workbook.add_format({"bg_color": "#ff0000"})
            worksheet.set_column(1, num_columns, 11)
            worksheet.conditional_format('G1:G1048576', {'type': '3_color_scale', 'min_type': 'num', 'min_value': 0,'mid_type': 'num','mid_value': 0.5,'max_type': 'num','max_value': 1})
            worksheet.conditional_format('I1:I1048576', {'type': 'cell', 'criteria': 'less than or equal to', 'value':  -0.65, 'format':   format_yellow})
            column_length = str(len(df))
            area_to_color = 'F1:' + "F" + column_length 
            worksheet.conditional_format(area_to_color, {'type': 'cell', 'criteria': 'equal to', 'value':  0, 'format':   format_red})
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)
            #Add Background colors
            bg_format1 = workbook.add_format({'bg_color': '#EEEEEE'}) # grey cell background color
            bg_format2 = workbook.add_format({'bg_color': '#FFFFFF'}) # white cell background color 
            for rn in range(len(df)):
                worksheet.set_row(rn, cell_format=(bg_format1 if rn%2==0 else bg_format2))
        #Close Writer
        writer.close()
        #Extract Output
        processed_data=output.getvalue()
        return processed_data


