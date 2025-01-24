import pytest
from unittest.mock import patch
from app.Mutations import MutationsModule
from app.VariantModel import Variant

# inputs = "inputs/"
# outputs = "outputs/"
# alignment_tool = "diamond"
# working_directory = os.getcwd()

# Run all tests with
# pytest test_5.py -v -rxs --color=auto --durations=0
# or
# pytest test_.py -v -rxs --color=auto --durations=0 -k "protein"

class TestMutations():
	tester_info = MutationsModule("/Users/karynmukiri/Desktop/rgi_total_resistome/tests/karyn_test/GCF_905400155.1_R-75386_assembly_genomic.fna.gz.temp.uncompressed.fsa.temp.contig.fsa.blastRes.xml", "contig", "GCF_905400155.1_R-75386_assembly_genomic.fna.gz", "/Users/karynmukiri/Desktop/rgi_total_resistome/tests/karyn_test")

	# @patch("app.Mutations.MutationsModule.xml_gen")
	def test_baseline(self):
		# print(MutationsModule.baseline(TestMutations.tester_info))

		"""
		old assert stuff
		"""
		pass
		# ## throwing dummy xml into our actual method (in mutations.py)
		# MutationsModule.xml_gen(xml_got)
		# ## our assert
		# mock_xml_gen.assert_called_with(xml_got)
		## just a print
	
		# ## initializing my MutationsModule class
		# mut_mod = MutationsModule(xml_got)

		# ## passing the dummy xml to the MAIN CLASS in our module (goes into __innit__)
		# tester_info = mut_mod(xml_got)

		# ## printing the output of our function of interest using our dummy xml as input
		# ## no patch required!
		# print(mut_mod.baseline(tester_info))
		
	def test_single_resistance_variant(self):
		Variant.run(TestMutations.tester_info)
		# print(MutationsModule.single_resistance_variant())
		# pass

	def test_VM(self):
		# print(Variant.run(TestMutations.tester_info))
		pass

	# def test_sum1(self, mocker):
	# 	# mocker.patch(__name__ + ".sum", return_value=5)
	# 	# assert sum(2, 3) == 5
	# 	pass