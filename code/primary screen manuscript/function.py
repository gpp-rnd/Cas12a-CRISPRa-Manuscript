import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.diagnostic import lilliefors
import pylab
def pvalue_normal_lilliefors(lfc_df_from_pdna, lfc_df_from_dropout, olfactory_gene_list):
    #inputs
    # lfc_df_from_pdna: LFC dataframe with pdna as reference; 
    #                   should have columns'Construct Barcode', 'Construct IDs', 
    #                   'untreated_repA', 'untreated_repB', 'avg_lfc_untreated'
    # lfc_df_from_dropout: LFC dataframe with untreated/dropout as reference; 
    #                   should have columns'Construct Barcode', 'Construct IDs', 
    #                   'Selumetinib_repA', 'Selumetinib_repB', 'avg_lfc_selumetinib'
    # olfactory_gene_list: a list of olfactory genes, they are genes with that starts with OR
    # Lilliefors_test was used to test normility because it does not specify the expected value and variance of the distribution
    
    #Dropout - pdna (Vailbility effect): all genes
    avg_lfc_untreated_test_statistics, avg_lfc_untreated_pval = lilliefors(lfc_df_from_pdna['avg_lfc_untreated'], 'norm')
    a_untreated_test_statistics, a_untreated_pval = lilliefors(lfc_df_from_pdna['untreated_repA'], 'norm')
    b_untreated_test_statistics, b_untreated_pval = lilliefors(lfc_df_from_pdna['untreated_repB'], 'norm')

    #Dropout - pdna (Vailbility effect) olfactory
    lfc_df_from_pDNA_olfactory = lfc_df_from_pdna[lfc_df_from_pdna['Construct IDs'].isin(olfactory_gene_list)]
    avg_lfc_untreated_olfactory_test_statistics, avg_lfc_untreated_olfactory_pval = lilliefors(lfc_df_from_pDNA_olfactory['avg_lfc_untreated'], 'norm')
    a_untreated_olfactory_test_statistics, a_untreated_olfactory_pval = lilliefors(lfc_df_from_pDNA_olfactory['untreated_repA'], 'norm')
    b_untreated_olfactory_test_statistics, b_untreated_olfactory_pval = lilliefors(lfc_df_from_pDNA_olfactory['untreated_repB'], 'norm')

    #Drug - Dropout: all genes
    avg_lfc_selumetinib_test_statistics, avg_lfc_selumetinib_pval = lilliefors(lfc_df_from_dropout['avg_lfc_selumetinib'], 'norm')
    a_selumetinib_test_statistics, a_selumetinib_pval = lilliefors(lfc_df_from_dropout['Selumetinib_repA'], 'norm')
    b_selumetinib_test_statistics, b_selumetinib_pval = lilliefors(lfc_df_from_dropout['Selumetinib_repB'], 'norm')

    #Drug - Dropout: olfactory
    lfc_df_from_dropout_olfactory = lfc_df_from_dropout[lfc_df_from_dropout['Construct IDs'].isin(olfactory_gene_list)]
    avg_lfc_selumetinib_olfactory_test_statistics, avg_lfc_selumetinib_olfactory_pval = lilliefors(lfc_df_from_dropout_olfactory['avg_lfc_selumetinib'], 'norm')
    a_selumetinib_olfactory_test_statistics, a_selumetinib_olfactory_pval = lilliefors(lfc_df_from_dropout_olfactory['Selumetinib_repA'], 'norm')
    b_selumetinib_olfactory_test_statistics, b_selumetinib_olfactory_pval = lilliefors(lfc_df_from_dropout_olfactory['Selumetinib_repB'], 'norm')

    data1 = {
        "condition": ['avg_lfc_untreated (all genes; lfc = untreated - pdna)',
                        'untreated_repA (all genes; lfc = untreated - pdna)',
                        'untreated_repB (all genes; lfc = untreated - pdna)',
                        'avg_lfc_untreated (olfactory genes; lfc = untreated - pdna)',
                        'untreated_repA (olfactory genes; lfc = untreated - pdna)',
                        'untreated_repB (olfactory genes; lfc = untreated - pdna)',
                        
                        'avg_lfc_selumetinib (all genes; lfc = selumetinib - untreated)',
                        'selumetinib_repA (all genes; lfc = selumetinib - untreated)',
                        'selumetinib_repB (all genes; lfc = selumetinib - untreated)',
                        'avg_lfc_selumetinib (olfactory genes; lfc = selumetinib - untreated)',
                        'selumetinib_repA (olfactory genes; lfc = selumetinib - untreated)',
                        'selumetinib_repB (olfactory genes; lfc = selumetinib - untreated)'],

        "pvalue": [avg_lfc_untreated_pval,
                  a_untreated_pval,
                  b_untreated_pval,
                  avg_lfc_untreated_olfactory_pval,
                  a_untreated_olfactory_pval,
                  b_untreated_olfactory_pval,
                  
                  avg_lfc_selumetinib_pval,
                  a_selumetinib_pval,
                  b_selumetinib_pval,
                  avg_lfc_selumetinib_olfactory_pval,
                  a_selumetinib_olfactory_pval,
                  b_selumetinib_olfactory_pval]}
    pvalue_df = pd.DataFrame(data1)
    return(pvalue_df)
#         "test_statistics": [avg_lfc_untreated_test_statistics,
#                            a_untreated_test_statistics,
#                            b_untreated_test_statistics,
#                            avg_lfc_untreated_olfactory_test_statistics,
#                            a_untreated_olfactory_test_statistics,
#                            b_untreated_olfactory_test_statistics,
                           
#                            avg_lfc_selumetinib_test_statistics,
#                            a_selumetinib_test_statistics,
#                            b_selumetinib_test_statistics,
#                            avg_lfc_selumetinib_olfactory_test_statistics,
#                            a_selumetinib_olfactory_test_statistics,
#                            b_selumetinib_olfactory_test_statistics],



def library_rename(library_df, library_name):
    #input: library_df: library dataframe that contains 'Construct Barcode', 'Construct IDs', 'Selumetinib_repA',
    #   'Selumetinib_repB', 'avg_lfc_selumetinib', 'z_scored_entire_dataset',
    #   'z_scored_olfactory_gene', 'olfactory_gene'
    #       library_name; str, name of the library
    # rename columns using the name of the library
    new_library_df = library_df.copy()
    for col in new_library_df.columns:
        if col in ['Construct Barcode', 'Selumetinib_repA','Selumetinib_repB', 'avg_lfc_selumetinib', 'z_scored_entire_dataset',
           'z_scored_olfactory_gene', 'olfactory_gene']:
            new_library_df = new_library_df.rename(columns={col: str(library_name +'_'+ col)})
    return(new_library_df)

def library_rename_viability(library_df, library_name):
    #input: library_df: library dataframe that contains 'Construct Barcode', 'Construct IDs', 'Selumetinib_repA',
    #   'Selumetinib_repB', 'avg_lfc_selumetinib', 'z_scored_entire_dataset',
    #   'z_scored_olfactory_gene', 'olfactory_gene'
    #       library_name; str, name of the library
    # rename columns using the name of the library
    new_library_df = library_df.copy()
    for col in new_library_df.columns:
        if col in ['Construct Barcode', 'untreated_repA','untreated_repB', 'avg_lfc_untreated', 
           'z_scored_olfactory_gene']:
            new_library_df = new_library_df.rename(columns={col: str(library_name +'_'+ col)})
    return(new_library_df)




    