from eye_tracking_analysis_complete import get_values_sj

df = get_values_sj("/Users/aryansmac/Documents/PROJECT/TEEP/TEEP\'25/LLD_CTRL_MAIN/Pupil-size-and-microsaccade-analysis-main/new_data/mix_data/derivatives/prec/sub-02", dv='Pupil_l', tw=(0, 2000), iv=['group'])
print(df)
