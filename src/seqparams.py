#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys
from multiseq import CytosineContext

#i path vanno dati come parametri espliciti
# vanno messi i valori di default


class SeqParams(object): #param wrapper
    def __init__(self):
        self.__avail_params = {
            "num_processes": 2,
#            "input_file": None,
#            "output_dir": None, "output_name": None,
#            "temp_dir": None,
            "seq_mode": "single_end",
            "lib_mode": "directional",
            "chromosomes": None,
            "fragment_length": None,
            "read_length": None,
            "p_meth_cg": None,
            "p_meth_chg": None,
            "p_meth_chh": None,
            "p_mutation": 0,
            "p_snp": 0,
            "coverage": 2.1
        }
        self.__values_params = dict()


    def params_summary(self):
        """ """
        print("******** I/O")
        print("- Input FASTA: {}".format(self.input_file))
        print("- Output directory: {}".format(self.output_dir))
        print("******** Sequencing options")
        print("- Chromosomes: {}".format(self["chromosomes"] if self["chromosomes"] is not None else "all"))
        print("- Sequencing: {}".format(self.seq_mode))
        print("- Library: {}".format(self.lib_mode))
        print("- Coverage: {}x".format(self.coverage))
        print("- Fragment size: {}".format(self.fragment_size))
        print("- Read length: {}".format(self.read_length))
        print("******** Probabilities")
        for ctx, p in self["p_meth"].items():
            print("- Methylation for {} context: {}".format(ctx, p))
        print("- SNP: {}".format(self.snp_rate))
        print("- Sequencing error: {}".format(self.error_rate))


    # def load_configuration(self, params):
    #     """ Load parameters from a configuration file """
    #
    #     self.load_params(params)
    #
    #     with open(params.config, "r") as f:
    #         for line in f:
    #             param, value = [elem.strip() for elem in line.split("=")]
    #
    #             if param in self.__avail_params:
    #                 self.__values_params[param] = value if value != "" else None
    #
    #     return self

    @property
    def input_file(self):
        return self.__values_params["input_file"]

    @property
    def output_dir(self):
        return self.__values_params["output_dir"]

    @property
    def output_name(self):
        return self.__values_params["output_name"]

    @property
    def temp_dir(self):
        return self.__values_params["temp_dir"]

    @property
    def seq_mode(self):
        return self.__values_params["seq_mode"]

    @property
    def lib_mode(self):
        return self.__values_params["lib_mode"]

    @property
    def num_processes(self):
        return self.__values_params["num_processes"]

    @property
    def fragment_size(self):
        return self.__values_params["fragment_length"]

    @property
    def read_length(self):
        return self.__values_params["read_length"]

    @property
    def error_rate(self):
        return self.__values_params["p_mutation"]

    @property
    def snp_rate(self):
        return self.__values_params["p_snp"]

    @property
    def coverage(self):
        return self.__values_params["coverage"]

    @property
    def p_meth(self):
        return {
            CytosineContext.CG: self.__values_params["p_meth_cg"],
            CytosineContext.CHG: self.__values_params["p_meth_chg"],
            CytosineContext.CHH: self.__values_params["p_meth_chh"]
        }


    def load_params(self, params):
        """ Load parameters from command line """

        self["input_file"] = params.fasta_file
        self["output_dir"] = params.output_path
        self["temp_dir"] = params.temp_dir
        self["chromosomes"] = params.chr
        self["seq_mode"] = params.seq_mode
        self["lib_mode"] = params.lib_mode
        self["fragment_length"] = params.fsize
        self["read_length"] = params.rlength
        self["num_processes"] = params.processes
        self["p_meth_cg"] = params.p_cg
        self["p_meth_chg"] = params.p_chg
        self["p_meth_chh"] = params.p_chh
        self["p_meth"] = self.p_meth
        self["p_mutation"] = params.error_rate
        self["p_snp"] = params.snp
        self["coverage"] = params.coverage

        self.__check_params()
        return self

    # def generate_config_file(self, filename):
    #     """ Generate an empty configuration file """
    #
    #     with open(filename, "w") as f:
    #         for param in self.__avail_params:
    #             f.write("{} = ...\n".format(param))
    #
    #     return self

    def __check_params(self):
        """ Program explodes if mandatory parameters are missing """ #ma con quale inglese

        #check mandatory parameter - input file
        if "input_file" not in self.__values_params:
            raise Exception("Input file missing")
        elif not os.path.exists(self.__values_params["input_file"]):
            raise Exception("Input file {} does not exist".format(self.__values_params["input_file"]))
        #check mandatory parameter - output stuff
        if "output_dir" not in self.__values_params:# or "output_name" not in self.__values_params:
            raise Exception("Output stuff is wrong")
        elif not os.path.exists(self.__values_params["output_dir"]):
            raise Exception("Output directory does not exist")
        #temporary directory
        if "temp_dir" not in self.__values_params:
            self.__values_params["temp_dir"] = "/tmp/"
        elif not os.path.exists(self.__values_params["temp_dir"]):
            raise Exception("temporary directory does not exist")
        #num processes
        if "num_processes" in self.__values_params:
            if not self.__cast_to_int("num_processes"):
                raise Exception("num_processes must be an integer value greater than zero")
        else:
            self.__values_params["num_processes"] = self.__avail_params["num_processes"]
        #sequencing type
        if "seq_mode" not in self.__values_params:
            self.__values_params["seq_mode"] = self.__avail_params["seq_mode"]
        #library mode
        if "lib_mode" not in self.__values_params:
            self.__values_params["lib_mode"] = self.__avail_params["lib_mode"]
        #fragment
        if "fragment_length" not in self.__values_params:
            raise Exception("fragment length value missing")
        elif not self.__cast_to_int("fragment_length"):
            raise Exception("fragment_length must be an integer value greater than zero")
        #read
        if "read_length" not in self.__values_params:
            raise Exception("read length value missing")
        elif not self.__cast_to_int("read_length"):
            raise Exception("read length must be an integer value greater than zero")
        #meth_probs
        if "p_meth" not in self.__values_params:
            raise Exception("Methylation probabilities are missing")
        #mutation prob
        if "p_mutation" not in self.__values_params:
            raise Exception("Mutation probability is missing")
        elif not self.__check_probability(self.error_rate):
            raise Exception("Probability values must be between 0 and 1")
        #snp prob
        if "p_snp" not in self.__values_params:
            raise Exception("Sequencing error probability is missing")
        elif not self.__check_probability(self.snp_rate):
            raise Exception("Probability values must be between 0 and 1")

    def __cast_to_int(self, attr):
        """ Casta attribute to int. Return a boolean value indicating the success of the operation"""
        try:
            x = self.__values_params[attr] = int(self.__values_params[attr])
        except ValueError:
            x = 0
        return x > 0

    def __check_probability(self, p):
        return p >= 0. and p <= 1.

    def __getitem__(self, key):
        return self.__values_params[key]

    def __setitem__(self, key, value):
        self.__values_params[key] = value
