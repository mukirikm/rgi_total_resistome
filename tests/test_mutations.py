import pytest
from unittest.mock import patch
from app.Mutations import MutationsModule

# inputs = "inputs/"
# outputs = "outputs/"
# alignment_tool = "diamond"
# working_directory = os.getcwd()

# Run all tests with
# pytest test_5.py -v -rxs --color=auto --durations=0
# or
# pytest test_.py -v -rxs --color=auto --durations=0 -k "protein"

class TestMutations():

	@patch("app.Mutations.MutationsModule.xml_gen")
	def test_xml_gen(self, mock_xml_gen):
		## dummy xml file
		xml_got = "/Volumes/juicy_thesis_secrets/scripts/rgi_total_resistome/tests/inputs/GCF_008988725.1_ASM898872v1_genomic.fna.gz.temp.uncompressed.fsa.temp.contig.fsa.blastRes.xml"

		## throwing dummy xml into our actual method (in mutations.py)
		MutationsModule.xml_gen(xml_got)
		## our assert
		mock_xml_gen.assert_called_with(xml_got)
		## just a print
		print(MutationsModule.xml_gen(xml_got))

	def test_fake(self):
		print("Works!")

	def test_sum1(self, mocker):
		mocker.patch(__name__ + ".sum", return_value=5)
		assert sum(2, 3) == 5