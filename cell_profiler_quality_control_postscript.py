"""
@author: Christian Gabriel

This script performs quality control on CellProfiler analyzed datsets.

"""
import os
import re
import pandas as pd
import numpy as np


# define folders


def analyze_dataset(cell_file,
                    main_colors,
                    tracking_marker,
                    output_file=None,
                    min_len = 48,
                    peak_threshold = 1.5,
                    size_jump_threshold = 0.2):
    
    """
    Parameters
    ---------
    cell_file: str 
        csv_file path, cellprofiler output 
    output_file: str (optional)
        path of the output file 
        creates dir if not present yet
        if not provided, "post_script_output.xlsx" will be stored in 
            input folder
    main_colors: list 
        color_name(s) which contain signals to be analyzed
    tracking_marker: str 
        color_name of the tracking channel used in CellProfiler
    min_len: int
        minimum length of approved time-series
    peak_threshold: float 
        threshold to define a peak in tracking channel for cell division
    size_jump_threshold: float 
        percentage of nucleus area to be flagged as "jump"
        indicative of cell division or tracking errors
    """
    
    
    ### Helper Functions
    def create_output_folder(output_file):
        output_folder=("/").join(output_file.split("/")[:-1])+"/"
        
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        return output_folder
    
    def load_data(cell_file):
        
        def erase_number(name):
            pattern = re.compile('[0-9]+')
            if re.match(pattern, name.split("_")[-1]):
                name = "_".join(name.split("_")[:-1])
            return name    
          
        #loading cell data
        # Standardize column labelling. Delete time series that were 
        #too short (decision in cell profiler pipeline)
        cell_data = pd.read_csv(cell_file)
        
        #erases arbitrary numbers in columns sometimes given by cell_profiler       
        cell_data=cell_data.rename(columns=erase_number)
        
        #Delete time series that were too short  
        cell_data=cell_data.dropna(subset = ["TrackObjects_Label"])
        
        return cell_data
    
    def extract_raw_time_series(cell_data, colors):
        
        signals_raw = {}
        colors_raw_extracted = []
        
        for color in colors:
            try:
                columns = "TrackObjects_Label"
                index="ImageNumber"
                values=f'Intensity_MeanIntensity_{color}_channel_corrected'
                signals_raw[color] = cell_data.pivot(index=index, 
                                                     columns=columns, 
                                                     values=values)
                colors_raw_extracted.append(color)
            except:
                continue                                    
        print(f'  extracted raw data from channels: {colors_raw_extracted}')
        
        return signals_raw
    
    def extract_bg_subtracted_time_series(cell_data, colors):
        
        signals_bg_subtracted = {}
        colors_bg_extracted = []
           
        for color in colors:
            try: 
                
                columns = "TrackObjects_Label"
                index="ImageNumber"
                values=f'Intensity_MeanIntensity_{color}_corr_BG'
                signals_bg_subtracted[color] = cell_data.pivot(index=index, 
                                                     columns=columns, 
                                                     values=values)
                colors_bg_extracted.append(color)
            except:
                continue
        print(f'  extracted BG subtracted data from channels: '
              f'{colors_bg_extracted}')                                             
        return signals_bg_subtracted        
    
    def extract_size_time_series(cell_data):
        return cell_data.pivot(index="ImageNumber", 
                                columns="TrackObjects_Label", 
                                values="AreaShape_Area", )
    
    def calculate_size_differences(object_sizes):
        absoute_size_differences = object_sizes.diff()
        relative_size_differences = absoute_size_differences/object_sizes
    
        return relative_size_differences
    
    def get_jumps(time_series, cutoff):
        #returns series, jumps above are marked with -1 or 1, all others are set to 0 
        
        def max_min_null(value, max_cutoff, min_cutoff=None):
        
            #evaluates value, if above max_cutoff returns 1, 
            #    if below min_cutoff returns -1, 
            #    else returns 0
            
            if not min_cutoff:
                min_cutoff= -1 * max_cutoff
            if value >= max_cutoff:
                return 1
            elif value <=min_cutoff:
                return -1
            else:
                return 0

        
        return time_series.apply(max_min_null, args =(cutoff,))

    def detect_cell_divisions(tracking_marker, signals_raw, 
                              size_jumps, peak_threshold):
        #Define cell divisions: Peak in RFP signal 
        #AND drop subsequent in cell size
        def detect_peaks_in_df(df, threshold):
            
            # detect peaks in iRFP signal, indicating cell division 
            def smooth(series):
                return series.rolling(7, center=True).mean()

            def pick_peaks_in_timeseries(timeseries, threshold):
                mean = timeseries.mean()
                std = timeseries.std()
                peaks = timeseries.apply(lambda x: 1 
                                         if x > mean + (threshold*std) 
                                         else 0)
                return peaks

            #smoothen time_series
            df_smooth = df.apply(smooth)
            #subtract smoothened time_series to expose peaks
            df_raw_minus_smooth = df - df_smooth
            #detects peaks which are > peak_threshold*std
            peaks_df = df_raw_minus_smooth.apply(pick_peaks_in_timeseries, 
                                                 args=(threshold,))
            return peaks_df
        
        cell_divisions = {}
        
        # detects peaks in the tracking marker intensity 
        # depending on peak_threshold
        tracking_marker_peaks = detect_peaks_in_df(signals_raw[tracking_marker]
                                                   , peak_threshold)
        
        #create a list with time-points of supposed cell_divisions per cell
        for cell in signals_raw[tracking_marker]:
            cell_divisions[cell]=[]
            #collect time-points of tracking marker peak per cell
            peak_times = list(tracking_marker_peaks[cell]
                              [tracking_marker_peaks[cell]==1].index)
            #if consequtive timepoints are marked peak, only keep the first one
            peak_times_single = [time for time in peak_times 
                                 if time-1 not in peak_times]
            
            #cell_division is defined as peak in tracking marker 
            #with decrease in object area at time of peak 
            #or within the folowing 2 timepoints
            for peak_time in peak_times_single:
                if (-1 in list(size_jumps[cell].loc[peak_time:peak_time+2])):
                    cell_divisions[cell].append(peak_time)
        return cell_divisions
    
    def get_non_division_size_jumps(size_jumps, cell_divisions):
        non_division_size_jumps={}
        for cell in size_jumps:
            non_division_size_jumps[cell]=[]
            #list all size jump times per cell
            size_jump_times = list(size_jumps[cell][size_jumps[cell]!=0].index)
            for size_jump_time in size_jump_times:
                #list only size jumps that are not occuring at division or 
                #within the folowing 2 timepoints
                if not size_jump_time in (cell_divisions[cell]+ 
                                           [x+1 for x in cell_divisions[cell]]+
                                          [x+2 for x in cell_divisions[cell]]):
                    non_division_size_jumps[cell].append(size_jump_time)
        return non_division_size_jumps
    
    def get_close_divisions(division_list):
        close_divisions = []
        for cell, div_time_list in division_list.items():
            if len(div_time_list) > 1:
               for position, div_time in enumerate(div_time_list[:-1]):
                   if ((div_time_list[position+1] - 
                        div_time_list[position])<15):
                       close_divisions.append(cell)
                       break
        return close_divisions
    
    def normalize_datasets_to_mean(dataset_dict):
        normalized_dataset_dict = {}
        for color, dataframe in dataset_dict.items():
            normalized_dataset_dict[color] = dataframe/np.mean(dataframe, 
                                                               axis=0)
        return normalized_dataset_dict    
    
    def filter_datasets(dataset_dict, min_len, 
                        non_division_size_jumps, 
                        close_divisions):
        
        def filter_cells(dataframe, non_division_size_jumps, 
                         close_divisions, min_len):
            filtered_df = dataframe.copy(deep=True)

            for cell, flag_list in non_division_size_jumps.items():
                cell_time_series = dataframe[cell].dropna()
                for flag in flag_list:

                    if (cell_time_series.loc[flag+1:].shape > 
                        cell_time_series.loc[:flag-1].shape):
                        cell_time_series = cell_time_series.loc[flag+1:]
                    else:
                        cell_time_series = cell_time_series.loc[:flag-1]
                if cell_time_series.shape[0] < min_len:
                    filtered_df.drop(cell, axis =1, inplace=True)    
                elif cell in close_divisions:
                    filtered_df.drop(cell, axis =1, inplace=True) 
                else:

                    filtered_df[cell] = cell_time_series
            return filtered_df

                
        filtered_dataset_dict = {}
        for color, dataframe in dataset_dict.items():
            filtered_dataset_dict[color] = filter_cells(dataframe, 
                                                        non_division_size_jumps, 
                                                        close_divisions, 
                                                        min_len)
        return filtered_dataset_dict
    
    def smoothen_out_divisions_dataset(dataset_dict, cell_divisions):
        def smoothen_out_divisions(dataframe, cell_divisions):
            df = dataframe.copy(deep=True)
            for cell, division_list in cell_divisions.items():
                try:
                    #extrapolating intensity at timepoint of division and 
                    #the one after linearily from surrounding values
                    for division_time in division_list:
                        
                        new_y_div_time = (df[cell].loc[division_time-1] +
                                           (df[cell].loc[division_time+2]- 
                                            df[cell].loc[division_time-1])
                                           *0.33)
                        new_y_div_time_plus_1 = (df[cell].loc[division_time-1]+
                                           (df[cell].loc[division_time+2]- 
                                            df[cell].loc[division_time-1])
                                           *0.66)
                        #set values to extrapolated values
                        df[cell].loc[division_time] = new_y_div_time
                        df[cell].loc[division_time+1] = new_y_div_time_plus_1
                except:
                    #occurs if division is detected at the edge of timeseries,
                    #since no values for extrapolation exist here. 
                    continue
            return df
        
                
        smoothened_dataset_dict = {}
        for color, dataframe in dataset_dict.items():
            #extrapolating intensity at timepoint of division and the one after
            #linearily from surrounding values
            smoothened_dataset_dict[color] = (
                smoothen_out_divisions(dataframe, 
                                       cell_divisions))
        return smoothened_dataset_dict
    
    def export_data_tables(cell_data, 
                          signals_bg_subtracted,
                          signals_bg_subtracted_rel, 
                          signals_cleaned, 
                          signals_cleaned_rel,
                          signals_cleaned_smooth,
                          signals_cleaned_rel_smooth,
                          all_cells,
                          approved_cells,
                          cell_divisions,
                          non_division_size_jumps,
                          min_len,
                          output_file):
        
        data_dict = {"raw_all_cells": signals_bg_subtracted,
                     "raw_acpt_cells": signals_cleaned,
                     "raw_acpt_cells_smoothDiv": signals_cleaned_smooth,
                     "norm_all_cells": signals_bg_subtracted_rel,
                     "norm_acpt_cells": signals_cleaned_rel,
                     "norm_acpt_cells_smoothDiv": signals_cleaned_rel_smooth,
                    }
        
        def export_xslx_overview(data_dict, output_file, *args, **kwargs):
            def create_overview(all_cells, approved_cells, signals_cleaned_smooth, min_len):
                overview = [["cells", len(all_cells)],
                            ["required conseq. images", min_len],
                            ["approved_cells",len(approved_cells)],
                           ]
                for color in signals_cleaned_smooth.keys():
                    overview.append([f'mean {color} signal',signals_cleaned_smooth[color].mean().mean()])
                return overview
            
            def create_flags(all_cells, cell_divisions, non_division_size_jumps, approved_cells):
                cell_divisions_and_flags = []
                for cell in all_cells:
                    cell_divisions_and_flags.append({"cell_number":cell,
                                    "divisions": cell_divisions[cell],
                                    "flags": non_division_size_jumps[cell],
                                    "approved": cell in approved_cells})
                return cell_divisions_and_flags
            
            overview = create_overview(all_cells, approved_cells, signals_cleaned_smooth, min_len)
            flags = create_flags(all_cells, cell_divisions, non_division_size_jumps, approved_cells)
            
            # Create a Pandas Excel writer using XlsxWriter as the engine.
            writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
                        
            pd.DataFrame(overview).to_excel(writer, sheet_name='overview')
            pd.DataFrame(flags).to_excel(writer, sheet_name='div_flags')
                        
            # Write each dataframe to a different worksheet.
            for data_name, data in data_dict.items():
                for color, dataframe in data.items():
                    dataframe.to_excel(writer, sheet_name=f'{color}_{data_name}')
     
            writer.save()
    
            
            
            
        def export_localization_data(cell_data, output_file):
            localization_data = ["Location_Center_X", "Location_Center_Y","ImageNumber",
                                 "TrackObjects_Label"]
            loc_df = cell_data[localization_data]
            loc_file_name = output_file[:-5]+"_localization_data.csv"
            loc_df.to_csv(loc_file_name)

            
        
        export_xslx_overview(data_dict, output_file, all_cells, approved_cells, signals_cleaned_smooth, 
                             min_len, cell_divisions, non_division_size_jumps)
        export_localization_data(cell_data, output_file)
    
    
    ### MAIN
        
    colors = main_colors + [tracking_marker]
    
    #create output folder
    
    if not output_file:
        input_folder = os.path.dirname(cell_file)
        output_file = input_folder+"/post_script_output.xlsx"
    else:
        create_output_folder(output_file)
        if not output_file.endswith(".xlsx"):
            output_file = output_file + ".xlsx"
               
    #load data files that were genereated by Cell Profiler
    #Standardize column labelling. Delete time series that 
    #were too short (decision in cell profiler pipeline)  
    cell_data = load_data(cell_file)
    
    #extract fluorescence data in dictionaries, keys are the colors
    signals_raw = extract_raw_time_series(cell_data, colors)
    signals_bg_subtracted = extract_bg_subtracted_time_series(cell_data, 
                                                              colors)
     
    #extract object size time-series
    object_sizes = extract_size_time_series(cell_data)   
    
    # From object size, calculate how the relative size of the cell 
    # has changed compared to previous time point    
    size_differences = calculate_size_differences(object_sizes)
    
    #detect size jumps over threshold
    size_jumps = size_differences.apply(get_jumps, 
                                        args=(size_jump_threshold,), 
                                        axis=0)
     
    #get a list of supposed cell divisions, defined as a peak in tracking
    #marker followed by decease in area
    cell_divisions = detect_cell_divisions(tracking_marker, signals_raw, 
                                           size_jumps, peak_threshold)
    
    #List size jumps that are not connected with divisions. 
    ##Most cases either misstracking, miss-segmentation, or edge-effects
    non_division_size_jumps = get_non_division_size_jumps(size_jumps, 
                                                          cell_divisions)

    #identify cells with too close devisions (<15h between divsions), 
    ##since probably something went wrong (miss-identified divisions, 
    #misstracking)          
    close_divisions = get_close_divisions(cell_divisions)
    
    #filter for high-confident tracks 
    #(no unexplicable size jumps, no close divisions),
    #with predifined minimum length
    signals_cleaned = filter_datasets(signals_bg_subtracted, min_len, 
                                      non_division_size_jumps, close_divisions)
    
    # fluorescence intensity of dividing/rounding cells are prone to artefacts
    # replace intensities at time of division by linear extrapolation of
    # surrunding values before and after division
    signals_cleaned_smooth = smoothen_out_divisions_dataset(signals_cleaned, 
                                                            cell_divisions)

    #calculating relative values
    signals_bg_subtracted_rel = (
        normalize_datasets_to_mean(signals_bg_subtracted)) 
    signals_cleaned_rel = (
        normalize_datasets_to_mean(signals_cleaned) )
    signals_cleaned_rel_smooth = (
        normalize_datasets_to_mean(signals_cleaned_smooth) )
    
    #create a list of all cells
    all_cells = [cell for cell in signals_bg_subtracted[tracking_marker]]
    #create a list of approved cells
    approved_cells = [cell for cell in signals_cleaned[tracking_marker]]
    
    print(f'{"{:.1f}".format(len(approved_cells)/len(all_cells)*100)} '
          '% of cells approved')
    
    #export data in one excel file    
    export_data_tables(cell_data, 
                        signals_bg_subtracted,
                        signals_bg_subtracted_rel, 
                        signals_cleaned, 
                        signals_cleaned_rel,
                        signals_cleaned_smooth,
                        signals_cleaned_rel_smooth,
                        all_cells,
                        approved_cells,
                        cell_divisions,
                        non_division_size_jumps,
                        min_len,
                        output_file)
    print(f'data stored at: {output_file}')         



#define primary colors which contain signals to be analyzed
main_colors = ["RFP", ]

#define tracking channel color
tracking_color = "iRFP"

input_folder =  "C:/Users/gabrielc/Desktop/Python/Githubscript/test/"       
input_filename  = "measurement_all_dataCells.csv"
output_filename = "output"


analyze_dataset(input_folder + input_filename, 
                main_colors,tracking_color,
                output_file = input_folder+"results/"+output_filename,
                )