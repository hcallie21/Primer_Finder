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
        has_dimer = "TAGCATTAGC"
        self.assertFalse(primer_eval.contains_dimer(no_dimer), "Should be False")
        self.assertTrue(primer_eval.contains_dimer(has_dimer), "should be true")
    #test contains hairpin
    def test_contains_hairpin(self):
        hairpin = "GGGTTTGGGTTT" 
        trinucleotide_loop = "GCGCAGC"
        self.assertTrue(primer_eval.contains_hairpin(hairpin), "Should be true")
        self.assertTrue(primer_eval.contains_hairpin(trinucleotide_loop), "Should be true")

            
    #test contains self complementary 
        #i think this is redundant
    #test melting temp
    def test_calculate_meltingtemp(self):
        empty = ""
        random = "ATGCCGAATGCATGCGTAC"
        random2 = "AATGGCCACGA"
        self.assertEquals(primer_eval.calculate_melting_temperature(empty), 0.0, "expected 0.0")
        self.assertEquals(primer_eval.calculate_melting_temperature(random), 67.5, "Expected 67.5")
        self.assertEquals(primer_eval.calculate_melting_temperature(random2),46.0, "Expected 46.0")
    #test check repeats 

    #test starts/ends w gc seq

    #test is suitable for primer

    #test find suitable primer