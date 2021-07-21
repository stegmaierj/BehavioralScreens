%%
 % CellCountingGUI.
 % Copyright (C) 2019 J. Stegmaier
 %
 % Licensed under the Apache License, Version 2.0 (the "License");
 % you may not use this file except in compliance with the License.
 % You may obtain a copy of the Liceense at
 % 
 %     http://www.apache.org/licenses/LICENSE-2.0
 % 
 % Unless required by applicable law or agreed to in writing, software
 % distributed under the License is distributed on an "AS IS" BASIS,
 % WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 % See the License for the specific language governing permissions and
 % limitations under the License.
 %
 % Please refer to the documentation for more information about the software
 % as well as for installation instructions.
 %
 % If you use this application for your work, please cite the repository and one
 % of the following publications:
 %
 % TBA
 %
 %%

helpText = {'Up Arrow: Increase contrast (gamma correction)', ...
            'Down Arrow: Decrease contrast (gamma correction)', ...
            'Left Arrow: Go to previous image', ...
            'Right Arrow: Go to next image', ...
            'A: Visualize the image with correct aspect ratio', ...
            'H: Show this help dialog', ...
            'O: Zoom out to the original view', ...
            'CTRL + Mouse Wheel: Zoom in/out in a Google-Maps like behavior', ...
            'Double-Click: Add a detection at the cursor position', ...
            'Right-Click: Remove detection closest to the cursor position', ...
            '', ...
            'Hint: In case key presses show no effect, left click once on the image and try hitting the button again. This only happens if the window looses the focus.'};

helpdlg(helpText);