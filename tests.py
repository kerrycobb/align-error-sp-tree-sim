#! /usr/bin/env python3

import sys
import os
import unittest
import tempfile
import dendropy
# import pyvolve
import random

import simulate 
import error


class TestBaseClass(unittest.TestCase):

    def setUp(self):
        self.scripts_dir_path = "scripts"
        self.biopy_path = "biopy-0.1.8"
        file_stream = tempfile.NamedTemporaryFile(
                prefix="tests",
                delete=False)
        self.tmp_path = file_stream.name
        file_stream.close()
        self.tmp_dir = tempfile.mkdtemp(prefix="tests")
        self.rng = random.Random()

    # def tearDown(self):
    #     os.remove(self.tmp_path)
    #     align_dir = os.path.join(self.tmp_dir, "alignments")
    #     for file_name in ("alignment-0.fa", "rates-0.txt", "info-0.txt"):
    #         path = os.path.join(align_dir, file_name)
    #         if os.path.exists(path):
    #             os.remove(path)
    #     if os.path.exists(align_dir):
    #         os.rmdir(align_dir)
    #     os.rmdir(self.tmp_dir)

    @classmethod
    def expected_yule_tree_height(cls, ntips, birth_rate):
        height = 0
        for i in range(2, ntips + 1):
            height += float(1) / (i * birth_rate)
        return height

    @classmethod
    def expected_yule_tree_length(cls, ntips, birth_rate):
        return float(ntips - 1) / birth_rate

    def run_simulate_species_tree(self, n_sp, birth_rate, death_rate=0.0,
            pop_size_shape=1.0,pop_size_scale=1.0):
        sp_tree = simulate.simulate_species_tree(
            n_sp=n_sp,
            birth_rate=birth_rate,
            death_rate=death_rate,
            pop_size_shape=pop_size_shape,
            pop_size_scale=pop_size_scale,
            rng=self.rng)
        return sp_tree

    def run_simulate_gene_tree(self, sp_tree, n_tips_per_sp):
        gtrees = simulate.simulate_gene_tree(
            sp_tree=sp_tree,
            n_tips_per_sp=n_tips_per_sp,
            rng=self.rng)
        return gtrees

    def run_simulate_alignment(self, path, gene_tree, locus_length):
        simulate.simulate_alignment(
            path=path,
            rng=self.rng,
            gene_tree=gene_tree,
            locus_length=locus_length)

    def run_drop_singletons(self, input, output, prob):
        error.drop_singletons(
            input=input,
            output=output,
            prop=prob, 
            rng=self.rng)


class TestSimulateSpeciesTree(TestBaseClass):           

    def sim_species_tree_testing(self, n_sp, birth_rate, death_rate=0.0,
            pop_size_shape=1.0, pop_size_scale=1.0, nreps=20000, age_places=2,
            length_places=2, pop_size_places=2):
        root_ages = []
        tree_lengths = []
        pop_sizes = []
        for i in range(nreps):
            sp_tree = self.run_simulate_species_tree(
                n_sp=n_sp,
                birth_rate=birth_rate,
                death_rate=death_rate,
                pop_size_shape=pop_size_shape,
                pop_size_scale=pop_size_scale)
            # sys.stdout.write("Tree: {0}\n".format(
            #         sp_tree.as_string(schema="newick",
            #                 annotations_as_nhx=True,
            #                 suppress_annotations=False)))
            root_ages.append(sp_tree.seed_node.age)
            tree_lengths.append(sp_tree.length())
            for node in sp_tree:
                pop_sizes.append(node.edge.pop_size)
        self.assertEqual(nreps, len(root_ages))
        self.assertEqual(nreps, len(tree_lengths))
        n_branches = (2 * n_sp) - 1
        self.assertEqual(nreps * n_branches, len(pop_sizes))
        mean_age = sum(root_ages) / len(root_ages)
        expected_age = self.expected_yule_tree_height(n_sp, birth_rate)
        mean_length = sum(tree_lengths) / len(tree_lengths)
        expected_length = self.expected_yule_tree_length(n_sp, birth_rate)
        mean_size = sum(pop_sizes) / len(pop_sizes)
        expected_size = pop_size_shape * pop_size_scale
        self.assertAlmostEqual(mean_age, expected_age, places=age_places)
        self.assertAlmostEqual(mean_length, expected_length, places=length_places)
        self.assertAlmostEqual(mean_size, expected_size, places=pop_size_places)

    def test_two_species(self):
        self.sim_species_tree_testing(
            n_sp=2,
            birth_rate=2.0,
            death_rate=0.0,
            pop_size_shape=10.0,
            pop_size_scale=0.1)

    def test_two_species(self):
        self.sim_species_tree_testing(
            n_sp=10,
            birth_rate=10.0,
            death_rate=0.0,
            pop_size_shape=2.0,
            pop_size_scale=1.0)


class TestSimulateGeneTrees(TestBaseClass):

    def test_two_species_one_tip(self):
        nreps = 50000
        sp_tree = self.run_simulate_species_tree(
            n_sp=2,
            birth_rate=5.0,
            death_rate=0.0,
            pop_size_shape=10.0,
            pop_size_scale=0.01)
        sys.stdout.write("Tree: {0}\n".format(
            sp_tree.as_string(
                schema="newick",
                annotations_as_nhx=True,
                suppress_annotations=False)))
        sp_tree_root_age = sp_tree.length() / 2.0
        sp_tree_root_size = sp_tree.seed_node.edge.pop_size
        sys.stdout.write("Species tree root age: {0}\n".format(sp_tree_root_age))
        sys.stdout.write("Species tree root size: {0}\n".format(sp_tree_root_size))
        sys.stdout.write("Species tree pop sizes:\n")
        for node in sp_tree:
            sys.stdout.write("\t{0}\n".format(node.edge.pop_size))
        root_ages = []
        for rep in range(0, nreps):
            gene_tree = self.run_simulate_gene_tree(
                sp_tree=sp_tree,
                n_tips_per_sp=1)
            gene_tree.calc_node_ages()
            root_ages.append(gene_tree.seed_node.age)
        self.assertEqual(nreps, len(root_ages))
        mean_age = sum(root_ages) / len(root_ages)
        # Expected root age of gene tree is:
        #   time*mu_rate + 2*Ne_root*mu_rate
        # if time and Ne are scaled by the mutation rate, this becomes:
        #   time + 2*Ne_root
        # But, Ne in dendropy is the number of gene copies, so it is:
        #   time + Ne_root
        expected_age = sp_tree_root_age + (1.0 * sp_tree_root_size)
        sys.stdout.write("Expected gene tree root age: {0}\n".format(expected_age))
        sys.stdout.write("Mean gene tree root age: {0}\n".format(mean_age))
        self.assertAlmostEqual(mean_age, expected_age, places=3)


class TestSimulateAlignments(TestBaseClass):

    def test_two_tip_tree(self):
        root_height = 0.01
        tree_str = "[&R] (s1_1:{0},s0_1:{0}):0.0;".format(root_height)
        gene_tree = dendropy.Tree.get(data=tree_str, schema="newick")
        gene_tree.label = "0"
        locus_length = 400000
        alignment_path = os.path.join(self.tmp_dir, "alignment-0.phy")
        self.run_simulate_alignment(
            path=alignment_path,
            gene_tree=gene_tree,
            locus_length=locus_length)
        tips = ['1']
        taxa = ['s0', 's1']
        matrix = dendropy.DnaCharacterMatrix.get(
            path=alignment_path,
            schema="phylip")
        self.assertEqual(len(matrix), 2)
        self.assertEqual(len(matrix[0]), locus_length)
        num_diffs = 0
        for i in range(locus_length):
            if matrix[0][i] != matrix[1][i]:
                num_diffs += 1
        diffs_per_site = num_diffs / float(locus_length)
        self.assertAlmostEqual(diffs_per_site, 2 * root_height, places=3)

# class TestDropSingletons(TestBaseClass):           

#     def drop_singleton_testing(self, nreps=20000, prob=0.5):
    

    # 



if __name__ == '__main__':
    unittest.main()
