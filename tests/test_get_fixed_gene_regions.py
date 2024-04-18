import unittest
import random
import numpy as np
import pandas as pd
import pandas.testing as pd_testing
import sys
from os import path
sys.path.append("../src")
import func_get_fixed_gene_regions as fix


class TestUtils(unittest.TestCase):
    # set up shared data to test with
    @classmethod
    def setUpClass(cls):
        # create a fake overlaps
        # T1 is overlapped by 2 both of which are start or end outside gene
        # T2 is overlapped by 2 both of which are start or end within gene
        overlaps = {"TranscriptID":["TULIP1","PARKO","TULIP1","PARKO", "PARKO"], 
            "chr": [1, 4, 1, 4, 4],
            "Gene_Start": [20,126458901,20,126458901, 126458901], 
            "Gene_Stop": [1200, 126475169, 1200, 126475169, 126475169], 
            "Length": [1180, 16268, 1180, 16268, 16268],
            "strand": ["+", "-", "+", "-", "-"],
            "chr2": [5, 5, 5, 5, 5],
            "BidID":["chr1-dreg", "chr4-dreg", "chr1-tfit", "chr4-tfit", "chr4-tfit,dreg"], 
            "Bid_Start": [11, 126460000, 1040, 126474500, 126475168], 
            "Bid_Stop": [45, 126459000, 1300, 126475000, 126476000], 
            "overlap": [34, 1800, 160, 5, 1], 
            "DIFF": [10, 20, 30, 40, 50]}
        cls.overlaps = pd.DataFrame(data=overlaps)
        new_regions_results = {"chr":[4, 4, 4, 4, 1, 1],  
                               "source":[".",".",".",".",".","."], 
                               "feature":["gene","exon","exon","exon", 
                                          "gene", "exon"], 
                               "start":[126458901, 126458901.0, 126459000.0, 126475000.0, 20, 45.0], 
                               "end":[126475169, 126460000.0, 126474500.0, 126475168.0, 1200, 1040.0], 
                               "score":[".",".",".",".",".","."], 
                               "strand":["-", "-", "-", "-", "+", "+"], 
                               "frame":[".",".",".",".",".","."], 
                               "attribute":['gene_id "PARKO"', 
                                           'gene_id "PARKO"; transcript_id "PARKO:1"', 
                                           'gene_id "PARKO"; transcript_id "PARKO:2"', 
                                           'gene_id "PARKO"; transcript_id "PARKO:3"', 
                                           'gene_id "TULIP"', 
                                           'gene_id "TULIP"; transcript_id "TULIP:1"']}
        cls.new_regions_results = pd.DataFrame(new_regions_results)
        ex_results = {"GeneID":["PARKO:1", "PARKO:2", "PARKO:3", "TULIP1:1"], 
              "Chr":["4", "4", "4", "1"], 
              "Start":[126458901.0, 126459000.0, 126475000.0, 45.0], 
              "End":[126460000.0, 126474500.0, 126475168.0, 1040.0], 
              "Strand":["-", "-", "-", "+"]}
        ex_results_swap = {"GeneID":["TULIP1:1", "PARKO:1", "PARKO:2", "PARKO:3"], 
              "Chr":["1", "4", "4", "4"], 
              "Start":[45.0, 126458901.0, 126459000.0, 126475000.0], 
              "End":[1040.0, 126460000.0, 126474500.0, 126475168.0], 
              "Strand":["+", "-", "-", "-"]}
        cls.overlaps_results = pd.DataFrame(data=ex_results)
        cls.overlaps_results_swap = pd.DataFrame(data=ex_results_swap)
        cls.gene_start = 5
        cls.gene_stop = 20
        cls.sorted_starts1 = [2, 19]
        cls.sorted_ends1 = [8, 40]
        cls.sorted_starts2 = [7, 19]
        cls.sorted_ends2 = [9, 40]
        cls.sorted_starts3 = [7, 15]
        cls.sorted_ends3 = [9, 17]
        cls.chrom = 2
        cls.transcript_name = "TT"
        cls.strand = "+"
    
    @classmethod
    def tearDownClass(cls):
        cls.overlaps = None
        cls.gene_start = None
        cls.gene_stop = None
        cls.sorted_starts1 = None
        cls.sorted_ends1 = None
        cls.sorted_starts2 = None
        cls.sorted_ends2 = None
        cls.sorted_starts3 = None
        cls.sorted_ends3 = None
        cls.chrom = None
        cls.transcript_name = None
        cls.strand = None
    
    # def test_get_non_overlap_regions(self):
    #     expected_results = [(8, 19)]
    #     self.assertEqual(expected_results, 
    #                      fix.get_non_overlap_regions(self.gene_start, 
    #                                                  self.gene_stop, 
    #                                                  self.sorted_starts1, 
    #                                                  self.sorted_ends1))
    #     expected_results = [(5, 7), (9, 19)]
    #     self.assertEqual(expected_results, 
    #                      fix.get_non_overlap_regions(self.gene_start, 
    #                                                  self.gene_stop, 
    #                                                  self.sorted_starts2, 
    #                                                  self.sorted_ends2))
    #     expected_results = [(5, 7), (9, 15), (17, 20)]
    #     self.assertEqual(expected_results, 
    #                      fix.get_non_overlap_regions(self.gene_start, 
    #                                                  self.gene_stop, 
    #                                                  self.sorted_starts3, 
    #                                                  self.sorted_ends3))
    # def test_get_non_overlap_regions_gtf(self):
    #     expected_results = pd.DataFrame({"chr":["2", "2"], "source":[".","."], 
    #                                      "feature":["gene","exon"],
    #                                      "start":[5,8], "end":[20,19], 
    #                                      "score":[".","."], "strand":["+","+"], "frame":[".","."], 
    #                                      "attribute":["gene_id TT", "gene_id TT; transcript_id TT:1"]})
    #     print(expected_results)
    #     print(fix.get_non_overlap_regions_gtf(self.transcript_name, 
    #                                                     self.chrom, 
    #                                                     self.strand, 
    #                                                     self.gene_start, 
    #                                                  self.gene_stop, 
    #                                                  self.sorted_starts1, 
    #                                                  self.sorted_ends1, [10, 20]))
    #     pd_testing.assert_frame_equal(expected_results, 
    #                      fix.get_non_overlap_regions_gtf(self.transcript_name, 
    #                                                     self.chrom, 
    #                                                     self.strand, 
    #                                                     self.gene_start, 
    #                                                  self.gene_stop, 
    #                                                  self.sorted_starts1, 
    #                                                  self.sorted_ends1, [10, 20]))
#         expected_results = pd.DataFrame({"chr":["2", "2", "2"], "source":[".",".", "."], 
#                                          "feature":["gene","exon", "exon"],
#                                          "start":[5,5,9], "end":[20,7,19], 
#                                          "score":[".",".","."], "strand":["+","+","+"], "frame":[".",".","."], 
#                                          "attribute":["gene_id TT", "gene_id TT; transcript_id TT:1", "gene_id TT; transcript_id TT:2"]})
#         pd_testing.assert_frame_equal(expected_results, 
#                          fix.get_non_overlap_regions_gtf(self.transcript_name, 
#                                                         self.chrom, 
#                                                         self.strand, 
#                                                         self.gene_start, 
#                                                      self.gene_stop, 
#                                                      self.sorted_starts2, 
#                                                      self.sorted_ends2, [10,20]))
#         expected_results = pd.DataFrame({"chr":["4","4", "4", "4"], "source":[".",".", ".","."], 
#                                          "feature":["gene","exon","exon", "exon"],
#                                          "start":[126458901,126458901,126459000,126475000], "end":[126475169,126460000,126474500,126475168], 
#                                          "score":[".",".",".","."], "strand":["-", "-", "-","-"], "frame":[".",".",".","."], 
#                                          "attribute":["gene_id PARKO", "gene_id PARKO; transcript_id PARKO:1", "gene_id PARKO; transcript_id PARKO:2", "gene_id PARKO; transcript_id PARKO:3"]})
        
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("PARKO", 
#                                       4, "-", 126458901, 126475169, 
#                                       sorted([126460000, 126474500, 126475168]), 
#                                       sorted([126459000, 126475000, 126476000]), [10,20,30]))
#         # if Bidirectional completely covers gene, keep original gene but label as bidirectional
#         expected_results = pd.DataFrame({"chr":["chr3", "chr3"], "source":[".","."], 
#                                          "feature":["gene","exon"],
#                                          "start":[72663296,72663296], "end":[72663499,72663962], 
#                                          "score":[".","."], "strand":["-","-"], "frame":[".","."], 
#                                          "attribute":["gene_id LOC", "gene_id LOC; transcript_id LOC:BID"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("LOC", 
#                                       "chr3", "-", 72663296, 72663499, 
#                                       [72662842], [72663962], [0]))
#         # if 1 Bidirectional and on left side 
#         expected_results = pd.DataFrame({"chr":["chr3", "chr3"], "source":[".","."], 
#                                          "feature":["gene","exon"],
#                                          "start":[8441,8548], "end":[10267,10267], 
#                                          "score":[".","."], "strand":["-","-"], "frame":[".","."], 
#                                          "attribute":["gene_id LEFT", "gene_id LEFT; transcript_id LEFT:1"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("LEFT", 
#                                       "chr3", "-", 8441, 10267, 
#                                       [8354], [8548], [50]))
#         # if 1 bidirectional and on right side
#         expected_results = pd.DataFrame({"chr":["chr3", "chr3"], "source":[".","."], 
#                                          "feature":["gene","exon"],
#                                          "start":[6459,6459], "end":[7083,6862], 
#                                          "score":[".","."], "strand":["-","-"], "frame":[".","."], 
#                                          "attribute":["gene_id RIGHT", "gene_id RIGHT; transcript_id RIGHT:1"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("RIGHT", 
#                                       "chr3", "-", 6459, 7083, 
#                                       [6862], [7144], [50]))
#         # if overlaps completely
#         expected_results = pd.DataFrame({"chr":["chr3", "chr3"], "source":[".","."], 
#                                          "feature":["gene","exon"],
#                                          "start":[7807,7807], "end":[7912,7912], 
#                                          "score":[".","."], "strand":["-","-"], "frame":[".","."], 
#                                          "attribute":["gene_id ENVELOP", "gene_id ENVELOP; transcript_id ENVELOP:BID"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("ENVELOP", 
#                                       "chr3", "-", 7807, 7912, 
#                                       [7807], [8185], [-19]))
#         # if overlaps completely
#         expected_results = pd.DataFrame({"chr":["chr3", "chr3"], "source":[".","."], 
#                                          "feature":["gene","exon"],
#                                          "start":[7807,7807], "end":[7912,7912], 
#                                          "score":[".","."], "strand":["-","-"], "frame":[".","."], 
#                                          "attribute":["gene_id ENVELOP", "gene_id ENVELOP; transcript_id ENVELOP:BID"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("ENVELOP", 
#                                       "chr3", "-", 7807, 7912, 
#                                       [7806], [7912], [1]))
#         #if overlaps completely
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("ENVELOP", 
#                                       "chr3", "-", 7807, 7912, 
#                                       [7807], [7912], [0]))
#         expected_results = pd.DataFrame({"chr":["chr3", "chr3"], "source":[".","."], 
#                                          "feature":["gene","exon"],
#                                          "start":[4005,4301], "end":[4388,4388], 
#                                          "score":[".","."], "strand":["-","-"], "frame":[".","."], 
#                                          "attribute":["gene_id TEST", "gene_id TEST; transcript_id TEST:1"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("TEST", 
#                                       "chr3", "-", 4005, 4388, 
#                                       [4005], [4301], [30]))
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("TEST", 
#                                       "chr3", "-", 4005, 4388, 
#                                       [4005], [4301], [30]))
#         expected_results = pd.DataFrame({"chr":["chr3", "chr3", "chr3"], "source":[".",".","."], 
#                                          "feature":["gene","exon","exon"],
#                                          "start":[6692,6692,7304], "end":[8019,6700,8019], 
#                                          "score":[".",".","."], "strand":["-","-","-"], "frame":[".",".","."], 
#                                          "attribute":["gene_id TEST", "gene_id TEST; transcript_id TEST:1", "gene_id TEST; transcript_id TEST:2"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("TEST", 
#                                       "chr3", "-", 6692, 8019, 
#                                       [6700], [7304], [30]))
#         expected_results = pd.DataFrame({"chr":["chr3","chr3"], "source":[".","."], 
#                                          "feature":["gene","exon"],
#                                          "start":[2255,2255], "end":[3463,3001], 
#                                          "score":[".","."], "strand":["-","-"], "frame":[".","."], 
#                                          "attribute":["gene_id CITED", "gene_id CITED; transcript_id CITED:1"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("CITED", 
#                                       "chr3", "-", 2255, 3463, 
#                                       [3001], [3463], [20]))
#         expected_results = pd.DataFrame({"chr":["chr3","chr3"], "source":[".","."], 
#                                          "feature":["gene","exon"],
#                                          "start":[2255,2255], "end":[3463,3001], 
#                                          "score":[".","."], "strand":["-","-"], "frame":[".","."], 
#                                          "attribute":["gene_id ORF", "gene_id ORF; transcript_id ORF:1"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_gtf("ORF", 
#                                       "chr3", "-", 5715, 5904, 
#                                       [5693], [6369], [20]))
    
#     def test_get_non_overlap_regions_df(self):
#         expected_results = pd.DataFrame({"GeneID":["TT:1"], "Chr":["2"], 
#                                        "Start":[8], "End":[19],  "Strand":["+"]})
#         pd_testing.assert_frame_equal(expected_results, 
#                          fix.get_non_overlap_regions_df(self.transcript_name, 
#                                                         self.chrom, 
#                                                         self.strand, 
#                                                         self.gene_start, 
#                                                      self.gene_stop, 
#                                                      self.sorted_starts1, 
#                                                      self.sorted_ends1, [10, 20]))
#         expected_results = pd.DataFrame({"GeneID":["TT:1", "TT:2"], "Chr":["2", "2"], 
#                                        "Start":[5, 9], "End":[7, 19],  "Strand":["+", "+"]})
#         pd_testing.assert_frame_equal(expected_results, 
#                          fix.get_non_overlap_regions_df(self.transcript_name, 
#                                                         self.chrom, 
#                                                         self.strand, 
#                                                         self.gene_start, 
#                                                      self.gene_stop, 
#                                                      self.sorted_starts2, 
#                                                      self.sorted_ends2, [10,20]))
        
#         expected_results = pd.DataFrame({"GeneID":["PARKO:1", "PARKO:2", "PARKO:3"], "Chr":["4","4","4"], "Start":[126458901, 126459000, 126475000], 
#                                          "End": [126460000, 126474500, 126475168], "Strand":["-", "-", "-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("PARKO", 
#                                       4, "-", 126458901, 126475169, 
#                                       sorted([126460000, 126474500, 126475168]), 
#                                       sorted([126459000, 126475000, 126476000]), [10,20,30]))
#         # if Bidirectional completely covers gene, keep original gene but label as bidirectional
#         expected_results = pd.DataFrame({"GeneID":["LOC:BID"], "Chr":["chr3"], "Start":[72663296], "End":[72663962], "Strand":["-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("LOC", 
#                                       "chr3", "-", 72663296, 72663499, 
#                                       [72662842], [72663962], [0]))
#         # if 1 Bidirectional and on left side 
#         expected_results = pd.DataFrame({"GeneID":["LEFT:1"], "Chr":["chr3"], "Start":[8548], "End":[10267], "Strand":["-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("LEFT", 
#                                       "chr3", "-", 8441, 10267, 
#                                       [8354], [8548], [50]))
#         # if 1 bidirectional and on right side
#         expected_results = pd.DataFrame({"GeneID":["RIGHT:1"], "Chr":["chr3"], "Start":[6459], "End":[6862], "Strand":["-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("RIGHT", 
#                                       "chr3", "-", 6459, 7083, 
#                                       [6862], [7144], [50]))
#         # if overlaps completely
#         expected_results = pd.DataFrame({"GeneID":["ENVELOP:BID"], "Chr":["chr3"], "Start":[7807], "End":[7912], "Strand":["-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("ENVELOP", 
#                                       "chr3", "-", 7807, 7912, 
#                                       [7807], [8185], [-19]))
#         # if overlaps completely
#         expected_results = pd.DataFrame({"GeneID":["ENVELOP:BID"], "Chr":["chr3"], "Start":[7807], "End":[7912], "Strand":["-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("ENVELOP", 
#                                       "chr3", "-", 7807, 7912, 
#                                       [7806], [7912], [1]))
#         #if overlaps completely
#         expected_results = pd.DataFrame({"GeneID":["ENVELOP:BID"], "Chr":["chr3"], "Start":[7807], "End":[7912], "Strand":["-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("ENVELOP", 
#                                       "chr3", "-", 7807, 7912, 
#                                       [7807], [7912], [0]))
#         expected_results = pd.DataFrame({"GeneID":["TEST:1"], "Chr":["chr3"], "Start":[4301], "End":[4388], "Strand":["-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("TEST", 
#                                       "chr3", "-", 4005, 4388, 
#                                       [4005], [4301], [30]))
#         expected_results = pd.DataFrame({"GeneID":["TEST:1"], "Chr":["chr3"], "Start":[4301], "End":[4388], "Strand":["-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("TEST", 
#                                       "chr3", "-", 4005, 4388, 
#                                       [4005], [4301], [30]))
#         expected_results = pd.DataFrame({"GeneID":["TEST:1","TEST:2"], "Chr":["chr3", "chr3"], 
#                                          "Start":[6692, 7304], "End":[6700, 8019], "Strand":["-", "-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("TEST", 
#                                       "chr3", "-", 6692, 8019, 
#                                       [6700], [7304], [30]))
#         expected_results = pd.DataFrame({"GeneID":["CITED:1"], "Chr":["chr3"], 
#                                          "Start":[2255], "End":[3001], "Strand":["-"]})
#         pd_testing.assert_frame_equal(expected_results, fix.get_non_overlap_regions_df("CITED", 
#                                       "chr3", "-", 2255, 3463, 
#                                       [3001], [3463], [20]))

#     def test_get_nonoverlap_SAF(self):
#         [new_regions_df, stats_df] = fix.get_nonoverlap_SAF(self.overlaps)
#         #print(stats_df)
#         stats_df.index = [0,1]
#         new_regions_df.index = [0,1,2,3]
#         stats_results = pd.DataFrame({"TranscriptID":["PARKO", "TULIP1"], 
#                                           "Bidirectionals":[["chr4-dreg", "chr4-tfit", "chr4-tfit,dreg"], ["chr1-dreg", "chr1-tfit"]], 
#                                           "Num_Bids":[3, 2], 
#                                           "BP_remain":[14462, 986], 
#                                           "Fraction_Overlap":[0.111, .164],
#                                           "Mean_Overlap":[602.0, 97.0], 
#                                           "Median_Overlap":[5.0, 97.0]})
#         stats_results_swap = pd.DataFrame({"TranscriptID":["TULIP1", "PARKO"], 
#                                           "Bidirectionals":[["chr1-dreg", "chr1-tfit"], ["chr4-dreg", "chr4-tfit", "chr4-tfit,dreg"]], 
#                                           "Num_Bids":[2, 3], 
#                                           "BP_remain":[986, 14462], 
#                                           "Fraction_Overlap":[.164, 0.111],
#                                           "Mean_Overlap":[97.0, 602.0], 
#                                           "Median_Overlap":[97.0, 5.0]})
#         if list(stats_df.TranscriptID) == ["PARKO", "TULIP1"]:
#             pd_testing.assert_frame_equal(self.overlaps_results, new_regions_df)
#             pd_testing.assert_frame_equal(stats_results, stats_df)
#         else:
#             pd_testing.assert_frame_equal(self.overlaps_results_swap, new_regions_df)
#             pd_testing.assert_frame_equal(stats_results_swap, stats_df)
        
        
#         # # Check works with TSSs
#         # TSS_list = ["chr5-dreg", "chr4-tfit,dreg"]
#         # [new_regions_df, stats_df] = fix.get_nonoverlap_SAF(self.overlaps, 
#         #                                                    set(TSS_list))
#         # stats_df.index = [0,1]
#         # new_regions_df.index = [0,1,2,3]
#         # stats_results = pd.DataFrame({"TranscriptID":["PARKO", "TULIP1"], 
#         #                                   "Bidirectionals":[["chr4-dreg", "chr4-tfit", "chr4-tfit,dreg"], ["chr1-dreg", "chr1-tfit"]], 
#         #                                   "Num_Bids":[3, 2], 
#         #                                   "BP_remain":[14462, 986], 
#         #                                   "Fraction_Overlap":[0.111, .164],
#         #                                   "Mean_Overlap":[602.0, 97.0], 
#         #                                   "Median_Overlap":[5.0, 97.0], 
#         #                              "Bidirectionals_nonTSS":[["chr4-dreg", "chr4-tfit"], ["chr1-dreg", "chr1-tfit"]], 
#         #                                   "Num_Bids_nonTSS":[2,2]})
#         # stats_results_swap = pd.DataFrame({"TranscriptID":["TULIP1", "PARKO"], 
#         #                                   "Bidirectionals":[["chr1-dreg", "chr1-tfit"], ["chr4-dreg", "chr4-tfit", "chr4-tfit,dreg"]], 
#         #                                   "Num_Bids":[2, 3], 
#         #                                   "BP_remain":[986, 14462], 
#         #                                   "Fraction_Overlap":[.164, 0.111],
#         #                                   "Mean_Overlap":[97.0, 602.0], 
#         #                                   "Median_Overlap":[97.0, 5.0], 
#         #                                   "Bidirectionals_nonTSS":[["chr1-dreg", "chr1-tfit"], ["chr4-dreg", "chr4-tfit"]], 
#         #                                   "Num_Bids_nonTSS":[2,2]})
#         # if list(stats_df.TranscriptID) == ["PARKO", "TULIP1"]:
#         #     pd_testing.assert_frame_equal(self.overlaps_results, new_regions_df)
#         #     pd_testing.assert_frame_equal(stats_results, stats_df)
#         # else:
#         #     pd_testing.assert_frame_equal(self.overlaps_results_swap, new_regions_df)
#         #     pd_testing.assert_frame_equal(stats_results_swap, stats_df)
#     def test_get_overlap_bid_stats(self):
#         overlaps = {"TranscriptID":["TULIP1","PARKO","TULIP1","PARKO", "PARKO", "THIRD"], 
#             "chr": [4, 4, 4, 4, 4, 4],
#             "Gene_Start": [47120,50000,47120,50000, 50000, 48000], 
#             "Gene_Stop": [48100, 51000, 48100, 51000, 51000, 50100], 
#             "Length": [1080, 1000, 1080, 1000, 1000, 2100],
#             "strand": ["+", "-", "+", "-", "-", "-"],
#             "chr2": [5, 5, 5, 5, 5, 5],
#             "BidID":["chr1-dreg", "chr1-dreg2", "chr1-tfit", "chr1-tfit2", "chr1-dreg", "chr1-dreg"], 
#             "Bid_Start": [48090, 50110, 47130, 50000, 48090, 48090], 
#             "Bid_Stop": [50010, 50450, 47500, 50090, 50010, 50010], 
#             "overlap": [34, 340, 370, 90, 10, 1920], 
#             "Bid_Length": [1920, 440, 370, 90, 1920, 1920], 
#             "Frac_Over_Gene": [0.031, 0.440, 0.343, 0.09, 0.01, 0.914], 
#             "Frac_Over_Bid": [0.018, 0.77, 1, 0.111, 0.005, 1]}
#         overlaps = pd.DataFrame(overlaps)
#         ex_results = {"BidID":["chr1-dreg", "chr1-dreg2", "chr1-tfit", "chr1-tfit2"],  
#                       "Transcripts":[["TULIP1", "PARKO", "THIRD"], ["PARKO"],
#                                     ["TULIP1"], ["PARKO"]], 
#                    "Num_Transcripts":[3, 1, 1, 1], 
#                       "Strand_Trans": [["+","-","-"], ["-"], ["+"], ["-"]],
#                "Fraction_Overlaps_bid":[[0.018, 0.005, 1], [0.77], [1], [0.111]], 
#                    "Fraction_Overlaps_genes":[[0.031, 0.01, 0.914], [0.440], [0.343], [0.09]], 
#                    "Tot_Fraction_Overlap_bid":[1.023, 0.77, 1, 0.111]}
#         ex_results = pd.DataFrame(ex_results)
#         stats_df = fix.get_overlap_bid_stats(overlaps)
#         stats_df = stats_df.sort_values("BidID")
#         stats_df.index = [0,1,2,3]
        
#         pd_testing.assert_frame_equal(ex_results, stats_df)
        
        
                         
if __name__ == "__main__":
    unittest.main()

                         
                         
        
