class Room {
    constructor(x, y, w, h, n) {
        this.x = x;
        this.y = y;
        this.h = h;
        this.w = w;
        this.c = `hsla(${~~(Math.random() * 360)},100%,50%,0.4)`;
        this.id = n;
    }
}
const rv = (n) => Math.floor(Math.random() * n),
    rir = (a, b) => Math.floor(Math.random() * (b - a)) + a;
new Vue({
    data: {
        // canv:document.querySelector('#canv'),
        rooms: [],
        rw: { s: 20, l: 100 },
        rh: { s: 30, l: 100 },
        rn: 5,
        minSep: 10,
        dm: false//"dungeon" mode: draw realistically
    },
    mounted() {
        const c = document.querySelector("#canv"),
            ctx = c.getContext("2d");
        this.canv = ctx;
        const srcs = [{ name: 'wall', url: `./img/lightrock.jpg` }, { name: 'floor', url: `./img/darkrock.jpg` }]
        // const wallIm = new Image(),
        //       floorIm = new Image(),;
        // wallIm.src = `https://img.glyphs.co/img?q=85&w=900&src=aHR0cHM6Ly9zMy5tZWRpYWxvb3QuY29tL2Jsb2ctaW1hZ2VzL1NSVC1JbWFnZS0wOS5qcGc/bXRpbWU9MjAxODEwMTUxNDIwNTY=`;
        // floorIm.src = `https://images-wixmp-ed30a86b8c4ca887773594c2.wixmp.com/f/8503b3f5-8649-470d-8d30-4f68035dd572/d8ysa65-dae49a94-b240-456e-8ebd-9640ee7603c8.png?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJ1cm46YXBwOjdlMGQxODg5ODIyNjQzNzNhNWYwZDQxNWVhMGQyNmUwIiwiaXNzIjoidXJuOmFwcDo3ZTBkMTg4OTgyMjY0MzczYTVmMGQ0MTVlYTBkMjZlMCIsIm9iaiI6W1t7InBhdGgiOiJcL2ZcLzg1MDNiM2Y1LTg2NDktNDcwZC04ZDMwLTRmNjgwMzVkZDU3MlwvZDh5c2E2NS1kYWU0OWE5NC1iMjQwLTQ1NmUtOGViZC05NjQwZWU3NjAzYzgucG5nIn1dXSwiYXVkIjpbInVybjpzZXJ2aWNlOmZpbGUuZG93bmxvYWQiXX0.LPwB-GTEPoFPQWCEFAZmbi5Oo-BbCjKuxPf8FRAgurI`;
        const proms = srcs.map(s => {
            return new Promise(resolve => {
                const im = new Image();
                im.addEventListener('load', () => {
                    resolve({ im: im, l: s.name })
                });
                im.src = s.url;
            })
        });
        Promise.all(proms).then(r => {
            console.log(r);
            r.forEach(q => {
                // document.querySelector('#main').appendChild(q.img)
                this[q.l] = this.canv.createPattern(q.im, 'repeat')
            })
            this.makeRooms();
        })

    },
    methods: {
        isCollided(rect1, rect2) {
            return (
                (rect1.x) < (rect2.x + rect2.w + this.minSep) &&
                (rect1.x + rect1.w + this.minSep) > (rect2.x) &&
                (rect1.y) < (rect2.y + rect2.h + this.minSep) &&
                (rect1.y + rect1.h + this.minSep) > (rect2.y)
            );
        },
        drawMap() {
            //blank canvas first
            console.log(this.wall, this.floor)
            this.canv.fillStyle = this.dm ? this.wall : "#222";
            this.canv.fillRect(0, 0, 500, 500);
            this.rooms.forEach((ra) => {
                this.canv.fillStyle = this.dm ? this.floor : ra.c;
                this.canv.fillRect(ra.x, ra.y, ra.w, ra.h);
            });
        },
        grc(rm) {
            return {id:rm.id,center:[rm.x + rm.w / 2, rm.y + rm.h / 2]};
        },
        makeRooms(keepOffMap) {
            //randomly place initial rooms
            this.rooms = new Array(this.rn)
                .fill(1)
                .map(
                    (q, i) =>
                        new Room(
                            rv(500),
                            rv(500),
                            rir(this.rw.s, this.rw.l),
                            rir(this.rh.s, this.rh.l),
                            i
                        )
                );
            //now progressively "shift" rooms until no overlaps
            /*
            For each room, we look at all rooms (filtering out the same room by ID)
            and see if it overlaps ("collides") with this room
            If so, we separate them by moving the leftmost one more left, 
            and the rightmost one more right, and same for topmost/bottommost
            We continue doing this until there are either no more loops,
            OR the maxMoves counter reaches zero (likely indicating an infinite loop!)
            */
            let maxMoves = 10000, //so we don't loop-lock
                hasCollision = true;
            while (hasCollision && maxMoves--) {
                hasCollision = false;
                this.rooms.forEach((ra) => {
                    this.rooms.forEach((rb) => {
                        if (ra.id == rb.id) return; //don't check room against itself
                        if (this.isCollided(ra, rb)) {
                            hasCollision = true;
                            const cents = [this.grc(ra), this.grc(rb)],
                                yDiff = (cents[1][1] - cents[0][1]) / ((cents[1][0] - cents[0][0]));
                            if (ra.x > rb.x) {
                                ra.x++;
                                rb.x--;
                            } else {
                                ra.x--;
                                rb.x++;
                            }

                            if (ra.y > rb.y) {
                                ra.y++;
                                rb.y--;
                            } else {
                                ra.y--;
                                rb.y++;
                            }
                        }
                    });
                });
            }
            //rooms now "placed". Optionally  remove out-of-view rooms
            if(!keepOffMap){
                let h = this.canv.canvas.height,
                    w = this.canv.canvas.width;
                this.rooms = this.rooms.filter(rm=>rm.x>0 && rm.y>0 && (rm.x+rm.w)<w && (rm.y+rm.h)<h)
            }
            //draw stuff!
            this.drawMap();
            const paths = this.getPaths();
            console.log('paths',paths)
            
            this.graph = new Graph(paths);
            //for each Path, add a vertex to our graph (duplicates obvs ignored)
            paths.map(q=>q.id).forEach(v=>this.graph.addVertex(v));
            //now get all connections and add an edge for each connection
            this.getConnections(paths.map(p=>p.id)).forEach(c=>{
                this.graph.addEdge(...c)
            })
            //convert the graph connections to object, w/ distance
            // this.graph.edgeObjs(paths);
            // this.graph.conList
            // console.log('CONNECTIONS',this.graph.conList,'MAX CON W',Math.max(...this.graph.conList.map(q=>q.weight)))
            console.log('GRAFF',this.graph,this.graph.conList, this.graph.mst(this.graph.conList))
            this.drawPaths(paths);
            // this.drawBranches(this.doMst());
        },
        getConnections(p){
            const cons = [];
            let tempCons = null;
            for(let i=0;i<p.length-1;i++){
                tempCons = [p[i],p[i+1]].sort();
                if(tempCons[0]==tempCons[1] || !!cons.find(cc=>cc[0]==tempCons[0] && cc[1]==tempCons[1])){
                    continue;
                }
                cons.push(tempCons)
            }
            return cons;
        },
        getPaths() {
            const centerPoints = this.rooms.map(this.grc),
                delInds = Delaunay.triangulate(centerPoints,'center');
            // console.log('centerPoints', centerPoints, 'tris', delInds);
            //note: delInds are INDICES of centerpoints, not centerpoints themselves!
            return delInds.map(q=>centerPoints[q]);
        },
        drawPaths(p) {
            const fp = p.shift();
            this.canv.strokeStyle = '#ddd';
            this.canv.beginPath()
            this.canv.moveTo(...fp.center);
            for(let i=0;i<p.length;i++){
                this.canv.lineTo(...p[i].center);
            }
            this.canv.stroke()
        },
        doMst(){
            const t = new Tree(this.rooms);
            console.log('tree',t)
            return t.spans;
        }
    }
}).$mount("#main");
