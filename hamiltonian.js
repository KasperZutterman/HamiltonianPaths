/*
Implementation of Hamiltonian Path algorithm due to 
Nathan Clisby, July 2012.
Modified by E. Emberly Sept 2014--!>
https://clisby.net/projects/hamiltonian_path/hamiltonian_path_v1.html
http://www.sfu.ca/~eemberly/phys847/assignments/hamiltonian_paths.html

Comments about the Markov chain used to generate paths
* using backbiting move described in Secondary structures in long
compact polymers, PHYSICAL REVIEW E 74, 051801 Ã�â€˜2006, by Richard
Oberdorf, Allison Ferguson, Jesper L. Jacobsen and Jan\'e Kondev
* algorithm is believed to be ergodic, but this has not been proved.
* current implementation is not the most efficient possible, O(N) for N
step walks, which could be improved with more sophisticated data
structure
* heuristic used for decision that equilibrium distribution is being
sampled from. This heuristic is quite conservative, but not certain.
* currently using default random number generator. This should be `good
enough' for generating typical walks, but shouldn't be replied upon for
serious numerical work.
*/
function backbite(n, path) {
    var i, j;
    var x, y;
    var dx, dy;
    var xedge, yedge;
    var iedge, add_edge;
    var success;
    /*
    Choose one of the two endpoints at random.
    Identify whether it is on the edge or not
    Decide which edge will be added (if any)
    Traverse walk, updating path
    */
    var itemp = Math.floor(Math.random() * 2);
    //Pre-compute n*n
    var nsq = n * n;
    if (itemp == 0) {
        //start at end: path[0]
        x = path[0][0];
        y = path[0][1];
        xedge = ((x == 0) || (x == n - 1));
        yedge = ((y == 0) || (y == n - 1));
        if (xedge && yedge) {
            /*
            corner
            1/3 acceptance probability
            */
            add_edge = Math.floor(Math.random() * 3) - 2;
        }
        else if (xedge || yedge) {
            /*
            edge
            2/3 acceptance probability
            */
            add_edge = Math.floor(Math.random() * 3) - 1;
        }
        else {
            //interior
            add_edge = Math.floor(Math.random() * 3);
        }
        success = (add_edge >= 0);
        iedge = 0;
        i = 3;
        while (iedge <= add_edge) {
            /*
            dx = x - path[i][0];
            dy = y - path[i][1];
            if (dx*dx+dy*dy == 1)
            */
            dx = Math.abs(x - path[i][0]);
            dy = Math.abs(y - path[i][1]);
            if (dx + dy == 1) {
                //we have found an empty edge
                if (iedge == add_edge) {
                    /*
                    This is the edge we wish to add.
                    reverse the walk from 0 to i-1
                    */
                    var jlim = (i - 1) / 2;
                    for (j = 0; j < jlim; j++) {
                        let temp = path[j];
                        path[j] = path[i - 1 - j];
                        path[i - 1 - j] = temp;
                    }
                }
                iedge++;
            }
            /*
            Can increment i by 2 due to bipartite nature of square
            lattice
            Even better: can increment by larger steps, but still
            ensure that we never miss a neighbour
            */
            i += Math.max(2, dx + dy - 1);
            /*
            i += 2;
            i++;
            */
        }
    }
    else {
        //start at end: path[nsq-1]
        x = path[nsq - 1][0];
        y = path[nsq - 1][1];
        xedge = ((x == 0) || (x == n - 1));
        yedge = ((y == 0) || (y == n - 1));
        if (xedge && yedge) {
            /*
            corner
            1/3 acceptance probability
            */
            add_edge = Math.floor(Math.random() * 3) - 2;
        }
        else if (xedge || yedge) {
            /*
            edge
            2/3 acceptance probability
            */
            add_edge = Math.floor(Math.random() * 3) - 1;
        }
        else {
            //interior
            add_edge = Math.floor(Math.random() * 3);
        }
        success = (add_edge >= 0);
        iedge = 0;
        i = nsq - 4;
        while (iedge <= add_edge) {

            /*
            dx = x - path[i][0];
            dy = y - path[i][1];
            if (dx*dx+dy*dy == 1)
            */
            dx = Math.abs(x - path[i][0]);
            dy = Math.abs(y - path[i][1]);
            if (dx + dy == 1) {
                //we have found an empty edge
                if (iedge == add_edge) {
                    /*
                    This is the edge we wish to add.
                    reverse the walk from i+1 to n*n-1
                    */
                    var jlim = (nsq - 1 - i - 1) / 2;
                    for (j = 0; j < jlim; j++) {
                        let temp = path[nsq - 1 - j];
                        path[nsq - 1 - j] = path[i + 1 + j];
                        path[i + 1 + j] = temp;
                    }
                }
                iedge++;
            }
            /*
            Can decrement i by 2 due to bipartite nature of square
            lattice
            Even better: can increment by larger steps, but still
            ensure that we never miss a neighbour
            */
            i -= Math.max(2, dx + dy - 1);
            /*
            i -= 2;
            i--;
            */
        }
    }
    return success;
}

function path_to_string(n, path) {

    var i, j, ii;
    // forward path
    ii = (path[0][0]) * n + path[0][1];
    var path_string = "" + ii + "";
    for (i = 1; i < n * n; i++) {
        ii = (path[i][0]) * n + path[i][1];
        path_string = path_string + " " + ii;
    }
    path_string += "\n";
    // reverse path
    ii = (path[n * n - 1][0]) * n + path[n * n - 1][1];
    path_string += "" + ii + "";
    for (i = n * n - 2; i >= 0; i--) {
        ii = (path[i][0]) * n + path[i][1];
        path_string = path_string + " " + ii;
    }
    path_string += "\n";
    return (path_string);
}

function generate_hamiltonian_path(n, q) {
    //initialize path
    var path = new Array(n * n);
    var i, j;
    var nsuccess, nattempts;
    for (i = 0; i < n; i++) {
        if (i % 2 == 0) {
            for (j = 0; j < n; j++) {
                path[i * n + j] = [i, j];
            }
        }
        else {
            for (j = 0; j < n; j++) {
                path[i * n + j] = [i, n - j - 1];
            }
        }
    }
    /*
    Now we attempt to apply backbite move repeatedly
    Our stopping criterion is that we want the random
    walk to have `covered' the whole grid.
    20*n*n successful moves is clearly not enough,
    by inspection on 100x100 grid.
    Take 10*n*n*n. By inspection this is enough, but slow.
    Relevant time for equilibrium is the cover time for an nxn grid.
    This is O(n^2 log^2 n)
    So ... could take const. * n^2 * log^2 n
    By inspection, this does a good job, and is asymptotically faster
    than previous proposol of O(n^3)
    */
    nsuccess = 0;
    nattempts = 0;
    /*
    Constant factor is a guess which is based on appearance - could
    experiment with making factor a bit smaller than 10.0, e.g. 5.0
    for faster run time, or maybe doubling it to 20.0 to ensure that
    the resulting path is truly random.
    For this reason, quality factor q introduced for 
    user to manipulate.
    */
    let nmoves = q * 10.0 * n * n * Math.pow(Math.log(2. + n), 2);
    /*
    nmoves = 10*n*n
    nmoves = 5*n*n*n
    */
    while (nattempts < nmoves) {
        let success = backbite(n, path);
        if (success) nsuccess++;
        nattempts++;
    }
    /*
    alert('Pr(success) = '+nsuccess/(0.0001+nattempts));

    Now we apply the same number of attempts.
    N.B.: if we just tested the number of successful moves then the
    result would be biased. (i.e. not truly 'random'), because it
    would be conditional on the last attempted move being successful

    for (i=0; i<nattempts; i++)
    {
        success = backbite(n,path);
    }
    */

    return path;
}

function generate_hamiltonian_circuit(n, q) {
    /*
    Generates circuits, but because we are subsampling circuits
    from the set of paths it is in fact not straightforward to
    sample uniformly at random from the set of circuits. Quite a subtle
    argument which I won't reproduce here.
    */
    var path = generate_hamiltonian_path(n, q);
    var nsq = n * n;
    var success;
    var min_dist = 1 + (n % 2);
    while (Math.abs(path[nsq - 1][0] - path[0][0])
        + Math.abs(path[nsq - 1][1] - path[0][1]) != min_dist) {
        success = backbite(n, path);
    }
    return path;
}

function draw_path(n, path) {
    var i;
    var x, y;
    var canvas = document.getElementById('path_canvas');
    if (canvas.getContext) {
        //Draw walk on the given canvas
        var ctx = canvas.getContext('2d');
        //clear canvas (this eliminates state information too)
        canvas.width = canvas.width;
        var w = canvas.width;
        var h = canvas.height;
        /*
        Calculate scale factors to convert to image location.
        Note that 2*pad is added to the denominator for padding.
        */
        var pad = 0.5;
        var sw = (w + 0.0) / (n + 2.0 * pad);
        var sh = (h + 0.0) / (n + 2.0 * pad);
        /*
        Choose line width. At the moment it is 
        0.2*(mean lattice spacing)
        */
        var lw = Math.floor(0.2 * 0.5 * (sw + sh));
        if (lw < 1) lw = 1;
        ctx.lineWidth = lw;
        //Draw path
        ctx.beginPath();
        ctx.moveTo((pad + path[0][0]) * sw, (pad + path[0][1]) * sh);
        var nsq = n * n;
        for (i = 1; i < nsq; i++) {
            ctx.lineTo((pad + path[i][0]) * sw, (pad + path[i][1]) * sh);
        }
        ctx.stroke();
        /*
        Draw endpoints as discs. The discs are chosen to be quite
        large so that they are visible for large n.
        */
        ctx.beginPath();
        x = (pad + path[0][0]) * sw;
        y = (pad + path[0][1]) * sh;
        ctx.arc(x, y, 0.4 * 0.5 * (sw + sh), 0.0, 2.0 * Math.PI, true);
        ctx.fillStyle = "rgb(255,0,0)";
        ctx.closePath();
        ctx.fill();
        ctx.beginPath();
        x = (pad + path[nsq - 1][0]) * sw;
        y = (pad + path[nsq - 1][1]) * sh;
        ctx.arc(x, y, 0.4 * 0.5 * (sw + sh), 0.0, 2.0 * Math.PI, true);
        ctx.fillStyle = "rgb(255,0,0)";
        ctx.closePath();
        ctx.fill();
    }
    return;
}

function draw_path_small(n, path, start_x, start_y, row_amount) {
    var i;
    var x, y;
    var canvas = document.getElementById('path_canvas');
    if (canvas.getContext) {
        //Draw walk on the given canvas
        var ctx = canvas.getContext('2d');
        //clear canvas (this eliminates state information too)
        //canvas.width = canvas.width;
        var w = canvas.width / row_amount;
        var h = canvas.width / row_amount;
        //var h = canvas.height / row_amount;
        /*
        Calculate scale factors to convert to image location.
        Note that 2*pad is added to the denominator for padding.
        */
        var pad = 0.5;
        var sw = (w + 0.0) / (n + 2.0 * pad);
        var sh = (h + 0.0) / (n + 2.0 * pad);
        /*
        Choose line width. At the moment it is 
        0.2*(mean lattice spacing)
        */
        var lw = Math.floor(0.2 * 0.5 * (sw + sh));
        if (lw < 1) lw = 1;
        ctx.lineWidth = lw;
        //Draw path
        //ctx.beginPath((pad + path[0][0] + start_x) * sw, (pad + path[0][1] + start_y) * sh);
        ctx.beginPath();
        ctx.moveTo((pad + path[0][0]) * sw + start_x * w, (pad + path[0][1]) * sh + start_y * w);
        //console.log('moveto', (pad + path[0][0] + start_x) * sw, (pad + path[0][1] + start_y) * sh);
        var nsq = n * n;
        for (i = 1; i < nsq; i++) {
            ctx.lineTo((pad + path[i][0]) * sw + start_x * w, (pad + path[i][1]) * sh + start_y * w);
            //console.log('lineto', (pad + path[i][0] + start_x) * sw, (pad + path[i][1] + start_y) * sh);
        }
        ctx.stroke();
        /*
        Draw endpoints as discs. The discs are chosen to be quite
        large so that they are visible for large n.
        */
       
        ctx.beginPath();
        x = (pad + path[0][0]) * sw  + start_x * w;
        y = (pad + path[0][1]) * sh  + start_y * w;
        ctx.arc(x, y, 0.4 * 0.5 * (sw + sh), 0.0, 2.0 * Math.PI, true);
        ctx.fillStyle = "rgb(255,0,0)";
        ctx.closePath();
        ctx.fill();
        ctx.beginPath();
        x = (pad + path[nsq - 1][0]) * sw + start_x * w;
        y = (pad + path[nsq - 1][1]) * sh + start_y * w;
        ctx.arc(x, y, 0.4 * 0.5 * (sw + sh), 0.0, 2.0 * Math.PI, true);
        ctx.fillStyle = "rgb(255,0,0)";
        ctx.closePath();
        ctx.fill();
        
    }
    return;
}

function draw_all_paths(n, paths, row_amount) {
    var canvas = document.getElementById('path_canvas');
    //var ctx = canvas.getContext('2d');
    canvas.width = window.innerWidth - 30;
    canvas.height = (canvas.width / row_amount + 1) * (paths.length / row_amount);
    console.log(canvas.width / row_amount * paths.length / row_amount);
    for(let i = 0; i < paths.length; i++) {
        draw_path_small(n, paths[i], i % row_amount, Math.floor(i / row_amount), row_amount);
    }
}

function refresh_path() {
    //var nstring = prompt("Grid size = ?","10");
    var Npath = parseInt(document.path_parameters.elements["Npath"].value); //Number of attempted paths
    var n = parseInt(document.path_parameters.elements["n"].value); //grid size
    var q = parseFloat(document.path_parameters.elements["q"].value); //quality factor
    var path;
    var paths = [];
    var ii;
    var path_string = "";
    var text_box = document.getElementById("path_text");
    for (let i = 1; i <= Npath; i++) {
        path = generate_hamiltonian_path(n, q);
        paths.push(path);
        ii = path[0][0] * n + path[0][1];
        path_string += path_to_string(n, path);
    }
    //text_box.value = path_string;
    //draw_path(n, path);
    var row_amount = parseInt(document.path_parameters.elements["r"].value); //Number of drawings per row
    if (Npath < row_amount) {
        row_amount = Npath;
    }
    draw_all_paths(n, paths, row_amount);
    return;
}

var button = document.getElementById("button");
button.addEventListener("click", refresh_path);