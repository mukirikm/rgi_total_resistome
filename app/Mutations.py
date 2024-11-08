# from app.Base import BaseModel
from app.settings import *
from app.VariantModel import *
from Bio.Blast import NCBIXML

class MutationsModule(object):
    """Class for mutation searches (SNVs, frameshifts, coSNPs, etc.)."""
	
    def __init__(self, xml_file):
        self.xml_file = xml_file
			
    def __repr__(self):
        """Returns Mutation class full object."""         
        return "Mutation({}".format(self.__dict__)

    def xml_gen(self):
        with open(self.xml_file, 'r') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                perfect = {}
                strict = {}
                loose = {}
                for alignment in blast_record.alignments:
                    align_title = alignment.title
                    orf_info = blast_record.query.encode('ascii','replace')

                    c = 0
                    barc = 0
                    for eachc in orf_info:
                        if barc >= 6:
                            break
                        elif eachc == '|':
                            barc += 1
                            c += 1
                        else:
                            c += 1
                    orf_from = orf_info[c:]

                    model_type_id = self.extract_nth_bar(align_title, 0)
                    # logger.info("model_type_id: {} ".format(model_type_id))
                    space_pos = align_title.index(' ')
                    hit_id = align_title[0:space_pos]
                    hit_id = hit_id.encode('ascii','replace')
                    model_descrpt =align_title[align_title.index(' ')+1:]
                    underscore_in_MD = model_descrpt.index('_')
                    model_id = model_descrpt[0:underscore_in_MD]
                    seq_in_model = model_descrpt[underscore_in_MD+1: model_descrpt.index(' ')]
                    pass_value = self.extract_nth_bar(alignment.title, 1)
                    # logger.info("pass_value: {}".format(pass_value))
	
    def single_resistance_variant(self):
        """Searches for SNVs in sequences."""
        evalue_snp = self.extract_nth_bar(self.align_title, 2)
        snpl = []
        snp_dict_list = []
        # temp = ""
        evalue_snp_dec = evalue_snp
        snpl = evalue_snp_dec.split(',')

        for each_snp in snpl:
            snp_dict_list.append({"original": each_snp[0], "change": each_snp[-1], "position": int(each_snp[1:-1])})

        for eachs in snp_dict_list:
            pos = eachs["position"]
            ori = eachs["original"]
            chan = eachs["change"]

            if self.hsp.sbjct_start < pos and (self.hsp.sbjct_start + self.real_sbjct_length) > pos:
                orf_protein_sequence = ""

                # if predicted_genes_dict:
                # 	if orf_info.strip() in predicted_genes_dict.keys():
                # 		orf_protein_sequence = str(Seq(predicted_genes_dict[orf_info.decode()], generic_dna).translate(table=11)).strip("*")
                # 	else:
                # 		orf_protein_sequence = str(Seq(predicted_genes_dict[orf_info.decode()[:orf_info.decode().index(' # ')]], generic_dna).translate(table=11)).strip("*")

                if self.predicted_genes_dict_protein:
                    if self.orf_info.strip() in self.predicted_genes_dict_protein.keys():
                        orf_protein_sequence = self.predicted_genes_dict_protein[self.orf_info.decode()].strip("*")
                    else:
                        orf_protein_sequence = self.predicted_genes_dict_protein[self.orf_info.decode()[:self.orf_info.decode().index(' # ')]].strip("*")

                if self.submitted_proteins_dict:
                    orf_protein_sequence = str(self.submitted_proteins_dict[self.orf_info.decode().split(" ")[0]])

                # logger.info("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
                # 			self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
                # 			+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

                # Report ONLY if the SNPs are present
                qry = int(pos) - self.hsp.sbjct_start + self.find_num_dash(self.hsp.sbjct, (int(pos) - self.hsp.sbjct_start))
                sbj = int(pos) - self.hsp.sbjct_start + self.find_num_dash(self.hsp.sbjct, (int(pos) - self.hsp.sbjct_start))

                if self.hsp.query[qry] == chan:
                    query_snps = {}
                    # logger.debug("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
                    # 		self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
                    # 		+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

                    # get position of mutation in the query sequence
                    d = int(pos) - self.hsp.sbjct_start - self.find_num_dash(self.hsp.query, (int(pos) - self.hsp.sbjct_start))
                    # print("hsp.sbjct_start: ", hsp.sbjct_start)
                    self.query_snps = {"original": ori, "change": chan ,"position": d+1}
                    
                    return self.query_snps