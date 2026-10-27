# Creating a Python Program to Print "Hello World"

To create the Python program, I installed **Git, Python, and VS Code** on my computer.

I then:

1. Created a new repository called `Assignments`.

2. Created a folder inside it called `day02` and a Python program file using GitHub.

3. Cloned the entire `Assignments` repository onto my computer using:

```bash
git clone URL
```

4. Edited the Python program file locally because I forgot that Python requires parentheses for the `print()` function.

5. Committed the changes and used:

```bash
git push
```

to sync the changes to GitHub.

At this point, a file called `.DS_Store` was also created. I didn't know what it was at first, but my commit worked.

6. I realized that my Python program file was not inside the `day02` folder, so I moved it using:

```bash
mv ./assignments/Python_program1.py ./assignments/day02
```

7. I committed the changes again.

The commit message showed:

```text
Changes to be committed:
    modified:   .DS_Store
    renamed:    Python_program1.py -> day02/Python_program1.py
```

This helped me understand why the `.DS_Store` file appeared.

8. Finally, I ran the Python program using:

```bash
python Python_program1.py
```

It printed:

```text
Hi there
```

because "Hello World" is outdated :)

9. I asked chat to help with formatting this file and he mentioned I can create a .gitignore to ignore DS_store which is a mac file that doesnt need to be on my github apparently
