#PURPOSE: test all functions in primer_eval
import unittest
import primer_eval


class TestPrimerCases(unittest.TestCase): 
    #test calc gc content
    def test_GC_content(self): 
        no_gc = "ATTTATTA"
        only_gc = "GGGCCCCC"
        hlf_gc = "GTCAGTCA"

        self.assertEqual(primer_eval.calculate_gc_content(no_gc), 0, "Should be 0")
        self.assertEqual(primer_eval.calculate_gc_content(only_gc), 100, "Should be 100")
        self.assertEqual(primer_eval.calculate_gc_content(hlf_gc), 50, "Should be 50%")
    #test contains dimer, right now dimer length is set to 2
    def test_contains_dimer(self): 
        no_dimer = "ATGCGATCGAT"
        has_dimer = "ATTACTGCATGCAATGCAG"
        self.assertFalse(primer_eval.has_self_dimer(no_dimer), "Should be False")
        self.assertTrue(primer_eval.has_self_dimer(has_dimer), "should be true")
    # #test contains hairpin
    # def test_contains_hairpin(self):
    #     no_hairpin = "AAAAAAAAAAA"
    #     hairpin = "GGGTTTGGGTTT"
    #     trinucleotide_loop = "GCGCAGC"
    #     self.assertFalse(primer_eval.contains_hairpin(no_hairpin), "always claiming true")
    #     self.assertTrue(primer_eval.contains_hairpin(hairpin), "Should be true")
    #     self.assertTrue(primer_eval.contains_hairpin(trinucleotide_loop), "Should be true")

            
    #test contains self complementary 
        #i think this is redundant
    #test melting temp
    def test_calculate_meltingtemp(self):
        empty = "" 
        random = "ATGCCGAATGCATGCGTAC"
        random2 = "AATGGCCACGACC"
        self.assertEqual(primer_eval.calculate_melting_temperature(empty), 0.0, 2)
        self.assertAlmostEqual(primer_eval.calculate_melting_temperature(random), 51.089,2)
        self.assertAlmostEqual(primer_eval.calculate_melting_temperature(random2),38.407, 2)
    #test check repeats 
    def test_containsRepeats(self): 
        no_repeats = "ATGCAGACT" #under 3 
        has_rep_same = "ATAGCATTTTGACTA"
        has_rep = "GGCGGCGGCGGC"

        self.assertFalse(primer_eval.contains_repeats(no_repeats))
        self.assertTrue(primer_eval.contains_repeats(has_rep))
        self.assertTrue(primer_eval.contains_repeats(has_rep_same))

    #test starts/ends w gc seq
    def test_startsends_GC(self):
        no_gc = "ATTAAGCCCGAAAAA"
        only_start = "ATGAATCATTACGAAAAA"
        only_end = "AAATAGGACCATGATC"
        both = "GCATCGGATTGCATC"
        self.assertFalse(primer_eval.starts_ends_with_gc(no_gc), "Should be False")
        self.assertTrue(primer_eval.starts_ends_with_gc(only_start), "should be true")
        self.assertTrue(primer_eval.starts_ends_with_gc(only_end), "should be true")
        self.assertTrue(primer_eval.starts_ends_with_gc(both), "Should be true")
    #test is suitable for primer
    
    #test find suitable primer

if __name__ == "__main__": 
    unittest.main()