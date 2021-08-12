#!/usr/bin/env python3

import csv
import os
import subprocess
import sys
import time

verbose = False


class Task:
    """
    Task class
    """

    def __init__(self, executable):
        self.executable = executable
        self.inputfiles = []
        self.outputfiles = []
        self.arguments = []

    def add_inputs(self, *args):
        """
        Method to add input files to the task
        """
        for arg in args:
            self.inputfiles.append(arg)

    def add_outputs(self, *args):
        """
        Method to add output files to the task
        """
        for arg in args:
            self.outputfiles.append(arg)

    def add_args(self, *args):
        """
        Method to add command-line arguments to the task
        """
        for arg in args:
            self.arguments.append(arg)

    def __str__(self):
        """
        Method to print the task
        """
        print("  * Task for executable " + self.executable)
        print("    - command-line: " + self.executable + ' '.join(self.arguments))
        print("    - input files: ")
        for file in self.inputfiles:
            print("        - " + file)
        print("    - output files: ")
        for file in self.outputfiles:
            print("        - " + file)

    def run(self):
        """
        Method to run the task
        """
        global verbose

        sys.stderr.write(
            "Running a " + self.executable + " task with input files {" + ', '.join(self.inputfiles) + "} " +
            "and output files {" + ', '.join(self.outputfiles) + "}\n")

        cmd = self.executable + " " + ' '.join(self.arguments)

        if verbose:
            redirect = None
            sys.stderr.write('\tRunning sub command: ' + cmd + "\n")
        else:
            redirect = subprocess.DEVNULL

        start = time.time()
        os.chdir("./data/20130502/")
        if subprocess.call(cmd, shell=True, stderr=redirect, stdout=redirect) != 0:
            sys.stderr.write('\tCommand ' + cmd + ' failed!')
            sys.exit(1)
        end = time.time()
        os.chdir("../../")
        sys.stderr.write("  [executed in " + str("{:.2f}".format(end - start)) + " seconds]\n")


class Workflow:
    """
    Workflow class
    """

    def __init__(self, name):
        self.name = name
        self.tasks = []
        self.files_to_download = {}

    def add_tasks(self, *tasks):
        """
        Method to add a task to the workflow
        """
        for task in tasks:
            self.tasks.append(task)

    def __str__(self):
        """
        Method to print the workflow, if needed
        """
        print("Workflow " + self.name + " with " + str(len(self.tasks)) + " tasks")
        for task in self.tasks:
            print(task)

    def run(self):
        """
        Method to execute the workflow sequentially
        """
        start = time.time()
        sys.stderr.write("Running the workflow sequentially...\n")

        # Make a dictionary of all the tasks' output files, some of which
        # serve as input to other tasks. The dictionary key is the file name,
        # and the value is true if the file has been produced already, false
        # otherwise. Note that this dictionary doesn't store entries to
        # the workflow's input file
        all_output_files = {}
        for task in self.tasks:
            for f in task.outputfiles:
                all_output_files[f] = False

        # Loop until tasks remain (yes, this loop "destroys" the workflow)
        while len(self.tasks) > 0:

            # Pick a ready task
            ready_task = None
            for task in self.tasks:
                ready = True
                # Check if all the task's inputfiles are available
                for file in task.inputfiles:
                    # if the input file is the output of another task
                    # check whether it has already been produced or not
                    if file in all_output_files:
                        if not all_output_files[file]:
                            ready = False
                # If the task was ready, we're done
                if ready:
                    ready_task = task
                break

            # If we found a ready task, we run it
            if ready_task is not None:
                # Run the task
                ready_task.run()
                # Mark its output files as produced
                for file in ready_task.outputfiles:
                    all_output_files[file] = True
                self.tasks.remove(ready_task)
            else:
                # This should never happen
                sys.stderr.write("FATAL ERROR: No ready task found\n")
                sys.exit(1)

        end = time.time()
        sys.stderr.write("Workflow execution done in " + str("{:.2f}".format(end - start)) + " seconds.")


if __name__ == "__main__":
    wf = Workflow("1000Genome")

    # Population Files
    populations = ["ALL", "AFR", "AMR", "EAS", "EUR", "GBR", "SAS"]

    f = open("data.csv")
    datacsv = csv.reader(f)
    step = 1000
    c_nums = []
    individuals_files = []
    sifted_files = []
    sifted_tasks = []

    for row in datacsv:
        base_file = row[0]
        threshold = int(row[1])
        counter = 1
        output_files = []

        # Individuals tasks
        c_num = base_file[base_file.find("chr") + 3:]
        c_num = c_num[0:c_num.find('.')]
        c_nums.append(c_num)

        while counter < threshold:
            stop = counter + step
            out_name = "chr{}n-{}-{}.tar.gz".format(c_num, counter, stop)
            output_files.append(out_name)

            j_individuals = Task("individuals.py")
            j_individuals.add_inputs(base_file, "columns.txt")
            j_individuals.add_outputs(out_name)
            j_individuals.add_args(base_file, c_num, str(counter), str(stop), str(threshold))
            wf.add_tasks(j_individuals)

            counter = counter + step

        # merge job
        j_individuals_merge = Task("individuals_merge.py")
        j_individuals_merge.add_args(c_num)

        for f_chrn in output_files:
            j_individuals_merge.add_inputs(f_chrn)
            j_individuals_merge.add_args(f_chrn)

        f_chrn_merged = "chr{}n.tar.gz".format(c_num)
        individuals_files.append(f_chrn_merged)
        j_individuals_merge.add_outputs(f_chrn_merged)

        wf.add_tasks(j_individuals_merge)
        # individuals_merge_jobs.append(j_individuals_merge)

        # Sifting Job
        f_sifting = "./sifting/{}".format(row[2])
        f_sifted = "sifted.SIFT.chr{}.txt".format(c_num)
        sifted_files.append(f_sifted)

        j_sifting = Task("sifting.py")
        j_sifting.add_inputs(f_sifting)
        j_sifting.add_outputs(f_sifted)
        j_sifting.add_args(f_sifting, c_num)

        wf.add_tasks(j_sifting)
        sifted_tasks.append(j_sifting)

    # Analyses jobs
    for i in range(len(individuals_files)):
        for f_pop in populations:
            # Mutation Overlap Job
            j_mutation = Task("mutation_overlap.py")
            j_mutation.add_args('-c', c_nums[i], '-pop', f_pop)
            j_mutation.add_inputs(individuals_files[i], sifted_files[i], f_pop, "columns.txt")
            f_mut_out = "chr%s-%s.tar.gz".format(c_nums[i], f_pop)
            j_mutation.add_outputs(f_mut_out)
            wf.add_tasks(j_mutation)

    # Run the workflow sequentially
    wf.run()
