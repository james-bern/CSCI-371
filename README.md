## [🐄 docs (click me) 🐌](https://github.com/james-bern/CSCI-371/wiki)

## quick start guide
<details>
  <summary>🍏 on your macbook</summary>

#### build and run
- **open terminal:** `⌘ + Space`; type `Terminal`; press `Enter`
- **(optional but recommended) old-style fullscreen Terminal:** hold `option`; click green circle
- **navigate to Desktop:** `cd Desktop`
- **install Mac build tools if you don't have them already**
  - **(option a) you alread have Xcode and are done!**
  - **(option b -- recommended if you don't have Xcode) install standalone command line tools:**
    - go to this random link https://developer.apple.com/download/all/?q=xcode
    - sign in (make an account if necessary)
    - download the command line tools for your operating system (to find out which OS you have, click the apple logo at the top left of the screen -> about this Mac)
  - **(option c) install Xcode:** https://apps.apple.com/us/app/xcode/id497799835?mt=12
- **clone codebase:** `git clone https://github.com/james-bern/CSCI-371.git`
  - **_if prompted, install developer tools and then repeat this step_**
- **navigate into repo:** `cd CSCI-371`
- **build and run / debug:** `source mac_build_and_run.sh --help`
#### debug
- **debug with VS Code:**
  - **install and configure VS Code**
    - **follow steps for 'Installation':** https://code.visualstudio.com/docs/setup/mac
    - **follow steps for 'Launching from the command line':** https://code.visualstudio.com/docs/setup/mac#_launching-from-the-command-line
  - **_install CodeLLDB when prompted_**
  - **(optional but recommended) install C/C++ extension:** https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools
- **❗ _if you randomly start having problems, you may need to shell uninstall and then shell install again_ (see 'Launching from the command line')**
#### misc
- **(recommended) set refresh rate to 60 Hertz**
  - **open System Preferences:** `⌘ + Space`; type `System Preferences`; press `Enter`
  - **click Displays**
  - Refresh Rate -> 60 Hertz
- **make a file executable (runnable as `./foo.bar`):** `chmod +x foo.bar`
- **open finder:** `open .`
- **show hidden folders in finder:** `⌘ + Shift + .`
- **show all file extensions in finder:** Finder -> Preferences... -> Advanced -> Show all filename extensions
</details>
<hr>
<details>
  <summary>🤷‍♂️ on your windows 10 computer (note: windows 11 may have multiple issues--e.g. Desktop may be located inside of OneDrive--see me in person and we'll make it work)</summary>

#### build and run
- **install git if you don't have it**
  - **install standalone (use all default options):** https://git-scm.com/download/win
- **install C++ build tools if you don't have them (any year should work)**
  - **option a (recommended) -- install standalone (:warning: you will need to scroll down; be careful you install the right thing):** https://visualstudio.microsoft.com/downloads/?q=build+tools
  - **option b -- install Visual Studio Community Edition:** https://visualstudio.microsoft.com/downloads/
- **open a x64 Native Tools Command Prompt:** press `⊞ Win`; type `x64`; press `Enter`
- **navigate to Desktop:** `cd "C:\Users\YOUR_USER_NAME\Desktop"`
- **clone codebase:** `git clone https://github.com/james-bern/CSCI-371.git`
- **navigate into repo:** `cd CSCI-371`
- **build and run / debug:** `windows_build_and_run.bat --help`
- **(optional but recommended) create a shortcut for the x64 Native Tools Command Prompt that starts in repo**
  - **bring up x64 Native Tools Command Prompt in start menu:** press `⊞ Win`; type `x64`
  - **bring up shortcut in file explorer:** click Open file location
  - **create a Desktop shortcut:** holding control, click and drag highlighted shortcut to your Desktop
  - **customize shortcut:**
    - right click Desktop shortcut -> change Start in field to `"C:\Users\YOUR_USER_NAME\Desktop\CSCI-371"`
    - go to General tab -> change name (no-name field) to `CSCI-371` (or whatever you would like)

#### debug
- **option a (so much better) -- debug with RemedyBG (vroom vroom)**
  - purchase here: https://remedybg.itch.io/remedybg
  - place remedybg.exe in CSCI-371 directory
- **option b (maybe easier at first, but trust me it's not worth it) -- debug with Visual Studio Code**
  - **install VS Code if you don't have it:** https://code.visualstudio.com/download
#### misc
- **open file explorer:** start .
- **show hidden files:** press `⊞ Win`; type Show hidden files; File Explorer Options -> Show hidden files... (also a good time to make sure Hide extensions for known file types is unchecked)
</details>
<hr>
<details>
  <summary>🐧 on a linux machine in the Ward Lab (TBL 301; has NUC's) or Unix Lab (TCL 312; computers with blue dots have GPU's)</summary>
 
#### build and run
- **open a terminal:** `CTRL + Alt + T`
- **(optional; breaks `Ctrl + V` binding for Vim) change Terminal's paste shortcut from `Ctrl + Shift + V` to `Ctrl + V`:** Terminal -> Preferences -> Shortcuts -> Edit -> Paste -> `CTRL + V`
- **navigate to Desktop:** `cd Desktop`
- **clone codebase:** `git clone https://github.com/james-bern/CSCI-371.git`
- **navigate into repo:** `cd CSCI-371`
- **build and run / debug:** `source ubuntu_build_and_run.sh --help`
#### debug
- **_install CodeLLDB when prompted_**
- **(optional but recommended) install C/C++ extension:** https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools
#### misc
- **make a file executable (runnable as `./foo.bar`):** `chmod +x foo.bar`
- **open file explorer:** `nautilus .`
- **show hidden folders in file explorer:** `Ctrl + H`
- **panic in terminal:** `Ctrl + C`
- **unfreeze the terminal (frozen by Ctrl + S):** `Ctrl + Q`
</details>

## more details
<details>
  <summary>🐛 debugging</summary>

#### basic usage
- click to the left of a line number (red circle should appear) to toggle a breakpoint
#### troubleshooting
- **known issue -- pause button doesn't work at first in VS Code**
  - workaround: hit a manual breakpoint and resume, then pause should work
</details>

<details>
 <summary>ℹ️ about the codebase</summary>
  
#### design philosophy
the CSCI-371 codebase was designed to
- compile very fast
- run fast enough
- be simple
- be small
- not require an IDE
 
#### fun facts
- the codebase is using a [unity build](https://hero.handmade.network/forums/code-discussion/t/2651-unity_builds_in_practice) (smoosh everything into one file and compile it; no incremental building; no relation to the game engine)

#### credits
- how to debug an executable from the command line with vscode https://docs.omnetpp.org/articles/omnetpp-debugging-tips/
</details>

<details>
 <summary>🎉 C/C++ fun facts</summary>

```c
// // stack-allocated fixed-size arrays

// stack-allocated fixed-size array of 16 doubles filled with garbage
double a[16];

// same thing but filled with zeroes
double a[16] = {}; 

// these three lines are equivalent
double a[3] = { 0, 0, 0 };
double a[3] = {};
double a[] = { 0, 0, 0 };

// consider a function foo that takes a pointer to a double
void foo(double *bar);
// if we call foo on a fixed size array bar, then inside the body of foo
// the argument bar will no longer be a fixed size array
// it will have "decayed into a pointer"
```

```c
// // heap-allocated arrays

// b is a pointer to heap-allocated array of 1000 doubles filled with garbage
double *b = (double *) malloc(1000 * sizeof(double));

// same thing but filled with zeroes
double *b = (double *) calloc(1000, sizeof(double));
// ^ equivalent to calling malloc as above and then calling
// memset(b, 0, 1000 * sizeof(double));       
```

```c
// // structs

// this is a struct containing two doubles followed by two ints followed by a pointer to a Poe    
struct Poe {                                                                               
    double x, y;                                                                             
    int i, j;                                                                              
    Poe *twin;                                                                             
};                                                                                         

// stack-allocated struct of type Poe filled with garbage
Poe poe;

// same thing but filled with zeroes
Poe poe = {};
// for the example struct, this is equivalent to Poe poe = { 0.0, 0.0, 0, 0, NULL };          
// which we can also write as Poe poe = { 0, 0, 0, 0, 0 } without the compiler complaining    

// stack-allocated fixed-size array of 5 Poe's filled with zeroes
Poe poes[5] = {};

// pointer to a heap-allocated array of 9999 Poe's filled with garbage
Poe *poes = malloc(9999 * sizeof(Poe));

// you can partially intialize an array or struct (unspecified arguments filled with zeros)  
// e.g. Poe poe = { 3.14, -4.2, 1 }; is the same as Poe poe = { 3.14, -4.2, 1, 0, NULL };    

// the journey continues...                                     
Poe poe_array[4] = { { 1.3, 0, 5 }, {}, { .2, .4, 1, 3 } };
```

```c
// // white space

void purple (int *cow) {               
     *cow += 5;                        
}                                      
                                       
// is the same as                         
                                       
void purple(int*cow){*cow+=5;}         
                                       
// is the same as                         
                                       
void                                   
   purple      (int                    
* cow) {                               
  *                                    
cow += 5     ;              }          
```

```c
// // struct / class basics

// both the . and the -> let me access the fields of a struct/class
// .  is when we have the thing                                    
// -> is when we have a pointer to the thing                       
// rule of thumb: pass by pointer when we need to modify the thing 
// (e.g. sbuff_push_back)                                          

struct Teacher {
    double height;
    char *mothers_maiden_name;
    int social_security_number;
};

int main() {
    Teacher jim = {};
    jim.social_security_number = 123;
    return 0;
}
```
</details>

