"""
File containing export functions of the CNVizard.
Authors: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
Company: University Hospital Aachen
Email: jerkrause@ukaachen.de
"""

import pandas as pd
import xlsxwriter
from io import BytesIO


class CNVExporter:
    """
    Class used to export the filtered DataFrame as Excel tables.
    """

    def __init__(self):
        """
        Constructor of the class CNVExporter.
        """
        pass

    def save_tables_as_excel(self, total_df: pd.DataFrame, bintest_df: pd.DataFrame, hom_del_df: pd.DataFrame, 
                             total_candi_df: pd.DataFrame, bintest_candi_df: pd.DataFrame, consecutive_del_df: pd.DataFrame, 
                             consecutive_dup_df: pd.DataFrame) -> bytes:
        """
        Writes the preset DataFrames to an Excel file, each preset represents an Excel sheet.

        Args:
            total_df (pd.DataFrame): .cnr DataFrame.
            bintest_df (pd.DataFrame): bintest DataFrame.
            hom_del_df (pd.DataFrame): .cnr DataFrame filtered for homozygous deletions.
            total_candi_df (pd.DataFrame): .cnr DataFrame filtered for candigenes.
            bintest_candi_df (pd.DataFrame): bintest DataFrame filtered for candigenes.
            consecutive_del_df (pd.DataFrame): .cnr DataFrame filtered for consecutive deletions.
            consecutive_dup_df (pd.DataFrame): .cnr DataFrame filtered for consecutive duplications.

        Returns:
            bytes: Output data passed to the Streamlit download button object.
        """
        output = BytesIO()
        writer = pd.ExcelWriter(output, engine="xlsxwriter")
        sheets = {
            "total": total_df,
            "bintest": bintest_df,
            "hom_del": hom_del_df,
            "total_candi": total_candi_df,
            "bintest_candi": bintest_candi_df,
            "consecutive_del": consecutive_del_df,
            "consecutive_dup": consecutive_dup_df
        }

        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name, header=False, index=False, startrow=1)
            worksheet = writer.sheets[name]
            header_format = writer.book.add_format({
                "bold": True,
                "text_wrap": True,
                "valign": "top",
                "fg_color": "#D7E4BC",
                "border": 1
            })
            format_yellow = writer.book.add_format({"bg_color": "#FFFF00"})
            format_red = writer.book.add_format({"bg_color": "#ff0000"})
            worksheet.set_column(1, len(df.columns), 11)
            worksheet.conditional_format('G1:G1048576', {
                'type': '3_color_scale',
                'min_type': 'num',
                'min_value': 0,
                'mid_type': 'num',
                'mid_value': 0.5,
                'max_type': 'num',
                'max_value': 1
            })
            worksheet.conditional_format('I1:I1048576', {
                'type': 'cell',
                'criteria': 'less than or equal to',
                'value': -0.65,
                'format': format_yellow
            })
            column_length = str(len(df))
            area_to_color = 'F1:' + "F" + column_length
            worksheet.conditional_format(area_to_color, {
                'type': 'cell',
                'criteria': 'equal to',
                'value': 0,
                'format': format_red
            })
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)

            bg_format1 = writer.book.add_format({'bg_color': '#EEEEEE'})
            bg_format2 = writer.book.add_format({'bg_color': '#FFFFFF'})
            for rn in range(len(df)):
                worksheet.set_row(rn, cell_format=(bg_format1 if rn % 2 == 0 else bg_format2))

        writer.close()
        return output.getvalue()

    def save_filtered_table_as_excel(self, filtered_df: pd.DataFrame, filtered_name: str) -> bytes:
        """
        Writes the filtered DataFrame to an Excel file.

        Args:
            filtered_df (pd.DataFrame): Filtered .cnr DataFrame.
            filtered_name (str): Name for the Excel sheet.

        Returns:
            bytes: Output data passed to the Streamlit download button object.
        """
        output = BytesIO()
        writer = pd.ExcelWriter(output, engine="xlsxwriter")
        sheets = {filtered_name: filtered_df}

        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name, header=False, index=False, startrow=1)
            worksheet = writer.sheets[name]
            header_format = writer.book.add_format({
                "bold": True,
                "text_wrap": True,
                "valign": "top",
                "fg_color": "#D7E4BC",
                "border": 1
            })
            format_yellow = writer.book.add_format({"bg_color": "#FFFF00"})
            format_red = writer.book.add_format({"bg_color": "#ff0000"})
            worksheet.set_column(1, len(df.columns), 11)
            worksheet.conditional_format('G1:G1048576', {
                'type': '3_color_scale',
                'min_type': 'num',
                'min_value': 0,
                'mid_type': 'num',
                'mid_value': 0.5,
                'max_type': 'num',
                'max_value': 1
            })
            worksheet.conditional_format('I1:I1048576', {
                'type': 'cell',
                'criteria': 'less than or equal to',
                'value': -0.65,
                'format': format_yellow
            })
            column_length = str(len(df))
            area_to_color = 'F1:' + "F" + column_length
            worksheet.conditional_format(area_to_color, {
                'type': 'cell',
                'criteria': 'equal to',
                'value': 0,
                'format': format_red
            })
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)

            bg_format1 = writer.book.add_format({'bg_color': '#EEEEEE'})
            bg_format2 = writer.book.add_format({'bg_color': '#FFFFFF'})
            for rn in range(len(df)):
                worksheet.set_row(rn, cell_format=(bg_format1 if rn % 2 == 0 else bg_format2))

        writer.close()
        return output.getvalue()

    def save_tables_as_excel_tsv(self, filtered_tsv: pd.DataFrame) -> bytes:
        """
        Writes the filtered .tsv DataFrame to an Excel file.

        Args:
            filtered_tsv (pd.DataFrame): Filtered .tsv DataFrame.

        Returns:
            bytes: Output data passed to the Streamlit download button object.
        """
        output = BytesIO()
        writer = pd.ExcelWriter(output, engine="xlsxwriter")
        sheets = {"filtered": filtered_tsv}

        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name, header=False, index=False, startrow=1)
            worksheet = writer.sheets[name]
            header_format = writer.book.add_format({
                "bold": True,
                "text_wrap": True,
                "valign": "top",
                "fg_color": "#D7E4BC",
                "border": 1
            })
            worksheet.set_column(1, len(df.columns), 11)
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)

            bg_format1 = writer.book.add_format({'bg_color': '#EEEEEE'})
            bg_format2 = writer.book.add_format({'bg_color': '#FFFFFF'})
            for rn in range(len(df)):
                worksheet.set_row(rn, cell_format=(bg_format1 if rn % 2 == 0 else bg_format2))

        writer.close()
        return output.getvalue()
