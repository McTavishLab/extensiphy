---

title: Command Line Tutorial

author: Jasper Toscani Field

output: html_document

---
# Command Line Tutorial

This tutorial introduces the complete basics of the command line. Its not designed to be an in-depth exposition on how to use all facets of the command line. Rather, this tutorial will get you started and give you some important information that you can follow up with your own research.

1. When you open up a terminal window, you'll notice that you don't really see anything. DONT PANIC. This is just a new way of interacting with your computer. Lets look around our new environment. The `ls` command lets you see what files and folders (directories) are in the same folder (directory) as you currently are.

```bash
$ls
Documents
Downloads
any_other_folders
any_files
```

2. When you use the `ls` command, you'll see files and directories specific to your computer. you can also use the `-a` flag to show hidden files and directories.

```bash
$ls -a
.
..
.bash_profile
.bashrc
Documents
Downloads
any_other_folders
any_files
```

Flags are used in a couple different ways. A flag can be used to tell the program your using that you wish to use some extra functionality (an option), such as when we told the `ls` command that we also wished to see hidden files by using the `-a` flag. Flags can also indicate a particular variable to a program, such as a necessary file. The times when a flag is used for one or the other use-case will rely on the the program's designer to explain.

3. The `pwd` command will become very imporant because it tells you the path to your current location. By path, I mean the navigation path through the directory structure of your computer to get to your current location. If you open a terminal and type the `pwd` command, you should see something like:

```bash
$pwd
/home/your_user_name
```

The `/` (slashes) indicate the transition between one directory (folder) and another directory (folder). So `your_user_name` is a directory and it can be found in the `home` directory.

To illustrate the importance of the `pwd` command, we'll introduce another basic command: `cd`. The `cd` command moves you into another directory. So if I open a terminal, want to look at the files and directories around me and then move into one of those directories, I could use the following sequence of commands.

```bash
$ls
Documents
Downloads
any_other_folders
any_files

$pwd
/home/your_user_name

$cd Documents

$pwd
/home/your_user_name/Documents
```

So we started in the directory `your_user_name` and saw that a directory `Documents` was next to us. We then changed directory (`cd`) into the `Documents` directory. Our new location was confirmed with the `pwd` command. Awesome!

4. An important aspect of the command line that we'll need to use a lot in bioinformatics is pointing our computer to files or directories that aren't in the same directory we are. To do this, we'll use some of the hidden files we revealed earlier. Lets use the `ls` command to look for a specific file that we suspect is in a different directory.

```bash
$ls /home/your_user_name/.bash_profile

/home/your_user_name/.bash_profile
```

We were able to point to a file in a different directory by giving the computer the path to that file. This is probably one of the most important aspects of using command line. When using a bioinformatics program (like Extensiphy), you will need to point your computer to various files and directories using this slash-directory-slash structure.

Hopefully you can see how using this pathing information will relate to other command line programs that need to be passed the location of multiple files or directories and use multiple flags to describe the desired actions you wish that program to make.
