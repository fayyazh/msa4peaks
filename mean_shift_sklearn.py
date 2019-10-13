#!/usr/bin/env python

import argparse
import numpy as np
from sklearn.cluster import MeanShift
import multiprocessing
#import multiprocessing.pool
from pybedtools import BedTool

information = """
Calculates the mean_shift/peaks for the given dataset.

Input:
* bed6 file containing dataset
* bandwidth (optional) to calculate the cluster radius/range
* default bandwidth is 5
* maximum distance (optional) for bedtools-merge, default is 500
* number of cores/processes (optional)
* Output file name (optional)

Output:
* sorted bed6 file with new_centroids/peaks for the given data set.
 default name "output.bed".

Example:
- reads input.bed file  with bandwidth 3, distance 300, number of processes 4
- and writes data_output.bed with new peaks/centroids in sorted order.
- mean_shift.py input.bed -b 3 -d 300 -p 4 -o data_output.bed

"""


class Mean_Shift():

    """
    bandwidth to perform clustering
    >>> ms = Mean_Shift(3)
    >>> ms.bandwidth
    3
    """

    def __init__(self, bandwidth):
        self.bandwidth = bandwidth

    def group_chrm(self, bed_file):
        # ####################################################################
        # group chromosomes with respect to name, strand and interval start
        # ####################################################################

        """
        groups bed_file data(ndarray) w.r.t chromosome and strand for centroids
        >>> ms = Mean_Shift(3) # bandwidth:  3
        >>> bd = [['chrX', '1', '+', '1'], ['chrX', '2', '+', '1'],
        ... ['chrX', '3', '+', '1'], ['chrX', '4', '+', '1'],
        ... ['chrX', '5', '+', '1'], ['chrX', '9', '-', '9'],
        ... ['chrX', '10', '-', '9'], ['chrX', '11', '-', '9'],
        ... ['chrX', '12', '-', '9'], ['chrX', '13', '-', '9'],
        ... ['chrY', '20', '+', '20'], ['chrY', '21', '+', '20'],
        ... ['chrY', '22', '+', '20'], ['chrY', '23', '+', '20'],
        ... ['chrY', '24', '+', '20']]
        >>> sorted(ms.group_chrm(bd).items())
        [(('chrX', '+', '1'), ['1', '2', '3', '4', '5']),\
 (('chrX', '-', '9'), ['9', '10', '11', '12', '13']),\
 (('chrY', '+', '20'), ['20', '21', '22', '23', '24'])]
        """

        chrm_details = dict()

        """
        reading bed_file to make unique keys  and
        eventually to have independent/unique dataset to process
        independently for centroids computation.
        """
        for b in bed_file:
            """
            reads chromosome name, strand of input_file and start
            from merged_data respectively, in oder to make unique
            key for chromosome_pos_val i.e ('chrX', '+', '1').
            """
            key_name_strand_pos = (b[0], b[2], b[3])

            """
            chromosome_start_position value of an input bed data.
            """
            chrm_pos_val = b[1]

            """
            adding elements to the dictionary
            key_name_strand_pos is used as key for
            the dictionary chrm_details.
            """
            if key_name_strand_pos in chrm_details:

                """
                append to an existing key one by one
                {('chrX', '+', '1') : ['1', '2']}
                {('chrX', '+', '1') : ['1', '2', '3']}
                {('chrX', '+', '1') : ['1', '2', '3', '4']}
                {('chrX', '+', '1') : ['1', '2', '3', '4', '5']}
                """
                chrm_details[key_name_strand_pos].append(chrm_pos_val)
            else:
                """
                create a new array of values
                {('chrX', '+', '1') : ['1']}
                {('chrX', '+', '1'): ['1', '2', '3', '4', '5'],
                ('chrX', '-', '9'): ['9']}
                """
                chrm_details[key_name_strand_pos] = [chrm_pos_val]

        return chrm_details

    def centroids_chrm(self, chroms_start):
        # ####################################################################
        # Determines chromosomes centres with MeanShift
        # ####################################################################

        """
        applying mean shift algorithm to compute centres w.r.t bandwidth(3)
        >>> ms = Mean_Shift(3)
        >>> x = np.array([1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 20, 21, 22,\
         23, 24])
        >>> data = ms.centroids_chrm(x)
        >>> sorted(data, key=lambda tup: tup[0])
        [array([3.]), array([11.]), array([22.])]
        """

        ms = MeanShift(bandwidth=self.bandwidth, bin_seeding=False,
                       min_bin_freq=1, cluster_all=True)

        """
        chroms_start is an array [1, 2, 3, 4, 5]
        and converts it into numpy array [[1, 2, 3, 4, 5]]
        required to make an ndarray as input for the fit()
        """
        chroms_start = np.array([chroms_start])

        """
        Performs clustering on input(ndarray, shape (n_samples, n_features))
        and stores the centroids in cluster_centers_
        """
        centroids = ms.fit(chroms_start.transpose()).cluster_centers_

        # centroids = [[3] [22]] [[11]]
        return centroids

    def clear_output_file(self, output_file):
        # ####################################################################
        # delete contents of output_file already exists
        # ####################################################################

        # delete the contents of an output file if already exists.
        with open(output_file, 'w'):
            pass

    def output_bed(self, centroids, output, chromosome):
        # ####################################################################
        # writes output to a bed file format
        # ####################################################################
        """
        writes to a bed_file
        >>> Mean_Shift(3).output_bed('1', "test.bed",
        ... [("chrX", '+', '1')]) is None
        True
        """

        # bed file fields
        name = 'X'
        score = '255'

        with open(output, "a") as wr:
            # key = [('chrY', '+', '22')] centroid = [array([[38.5], [22.5]])]
            for key, centroid in zip(chromosome, centroids):
                # chromosome =  'chrY', strand = '+', _ = '22'
                chromosome, strand, _ = key

                for i in range(len(centroid)):
                    chrm_start = int(centroid[i])
                    chrm_end = chrm_start + 1

                    bed_data = ("%s\t%s\t%s\t%s\t%s\t%s" %
                                (chromosome, chrm_start,
                                 chrm_end, name,
                                 score, strand))
                    wr.write(bed_data + '\n')

        # sorting and saving as bedfile
        BedTool(output).sort().saveas(output)

    def call_centroid_and_output(self, chrm_details, output_file, core):
        # ####################################################################
        # calls to centroids_chrm method to determine centroids
        # calls to output_bed method to write to an output file
        # ####################################################################
        """
        determined centroid(s) and chromosomes information is written to file
        >>> import collections
        >>> chrm =  {('chrX', '+', '1'): ['1', '2', '3', '4', '5'],
        ... ('chrX', '-', '9'): ['9', '10', '11', '12', '13'],
        ... ('chrY', '+', '20'): ['20', '21', '22', '23', '24']}
        >>> dt = collections.OrderedDict(sorted(chrm.items()))
        >>> Mean_Shift(3).call_centroid_and_output(dt, "test.bed", 3)
        number of estimated clusters in ('chrX', '+', '1') are 1
        number of estimated clusters in ('chrX', '-', '9') are 1
        number of estimated clusters in ('chrY', '+', '20') are 1
        """

        chrom_keys = []
        chrom_values = []

        # making list of keys and values of chrm_details.
        for k in chrm_details:
            """
            k = ('chrX', '+', '1')
            making a list of keys
            [('chrX', '+', '1'), ('chrY', '+', '20')]
            """
            chrom_keys.append(k)

            """
            converts retrieved chromosomes start positions list to int
            ['1', '2', '3', '4', '5'] -> [1, 2, 3, 4, 5]
            and making a list for Pool as a parameter
            [[1, 2, 3, 4, 5], [20, 21, 22, 23, 24]]
            """
            chrom_values.append(list(map(int, chrm_details[k])))

        # chrom_values = [[1, 2, 3, 4, 5], [20, 21, 22, 23, 24]]
        pool = multiprocessing.Pool(processes=core)
        centroids = pool.map(self.centroids_chrm, chrom_values)
        pool.close()
        pool.join()
        pool.terminate()
        for i in range(len(centroids)):
            print("number of estimated clusters in %s are %s" %
                  (chrom_keys[i], len(centroids[i])))

        self.output_bed(centroids, output_file, chrom_keys)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=information,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument("data_file", type=str, help='Enter the bed file name')
    parser.add_argument("-b", "--bandwidth", type=int, default=5,
                        help='Enter the bandwidth')
    parser.add_argument("-d", "--distance", type=int, default=500,
                        help='Enter the maximum distance for bedtools-merge')
    parser.add_argument("-p", "--processes", type=int, default=2,
                        help='Enter the cores to use')
    parser.add_argument("-o", "--output_file", default='output.bed',
                        help='Enter the output bed file name')

    args = parser.parse_args()

    # file to save output
    output_file = args.output_file

    """
    Maximum distance between features to be merged.
    default distance is 500
    """
    distance = args.distance

    """
    no. of cores for multiprocessing
    default cores are 2
    """
    core = args.processes

    # default bandwidth is 5
    ms = Mean_Shift(args.bandwidth)
    print("Bandwidth : ", args.bandwidth)

    # delete the contents of output file if already exist
    ms.clear_output_file(output_file)

    # segregating bed file for independent data features
    input_file = BedTool(args.data_file)

    """
    s=True for only merge features that are the same strand.
    c=(2,6) gives max start position of chromosomes being merged
    and strand due to operation max and first on respective columns
    max and first are retrieved to make bed format for intersect.
    """
    merged_data = BedTool(input_file.sort().
                          merge(s=True, c=(2, 6), d=distance,
                                o=('max', 'first')))

    """
    s=True, reports hits in 'b' that overlap 'a' on the same strand.
    loj=True (Left outer join). Report features in 'a' with
    and without overlaps.
    'a' is input file in bed format
    'b' merged file in bed format
    each feature in 'a' is compared to 'b' in search of overlaps
    intersected data contains 12 columns, 6 of input file
    and 6 merged_data it overlaps with.
    """
    intersected_data = merged_data.\
        intersect(s=True, loj=True, a=input_file, b=merged_data)

    """
    intersected_data = chrX  1  5  X  255  +  chrX  1  6  +  5  +
    retrieves chromosome name, start and strand from input bed file
    and chromosome start from database/merged_data respectively.
    """
    bed_data = []
    for data in intersected_data:
        bed_col = str(data).strip().split()
        bed_data.append([bed_col[0], bed_col[1], bed_col[5], bed_col[7]])

    # groups data based on chromosome, strand; and start of overlap
    print(bed_data)
    pool = multiprocessing.Pool(processes=core)
    chrm_details = pool.map(ms.group_chrm, [bed_data])
    pool.close()
    pool.join()
    pool.terminate()

    # calls method for centroids computation and writes output to a file
    ms.call_centroid_and_output(chrm_details[0], output_file, core)
