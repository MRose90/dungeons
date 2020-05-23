// class Vertex {
//     constructor(p) {
//         this.x = p.x+(p.w/2);
//         this.y = p.y+(p.h/2);
//         this.to = null;
//         this.id = p.id;
//     }
//     dist(b) {
//         //calc euclidian distance
//         return Math.sqrt(Math.pow(Math.abs(b.x - this.x),2) + Math.pow(Math.abs(b.y - this.y),2))
//     }
// }
// class Tree {
//     constructor(verts) {
//         this.verts = verts.map(q=>new Vertex(q));
//     }
//     get spans() {
//         let currVert = 0,//our vertex "counter"
//             vertsDone = 0,
//             loopHalter = 10000;//for now, so we dont do infinite loops
//         // while (vertsDone < this.verts.length && loopHalter--) {
//         // }
//         for(let i=0;i<this.verts.length;i++){
//             this.verts[i].to = this.verts.filter(q=>q.id!==this.verts[i].id).sort((a,b)=>{
//                 return this.verts[i].dist(a) - this.verts[i].dist(b);
//             })[0]
//         }
//         return this.verts;
//     }
// }

class Graph {
    // https://hackernoon.com/the-javascript-developers-guide-to-graphs-and-detecting-cycles-in-them-96f4f619d563
    constructor(p) {
        this.paths = p;
        this.adjList = {};
    }
    addVertex(v) {
        this.adjList[v] = []
    }
    dist(a, b) {
        //calc euclidian distance
        return Math.sqrt(Math.pow(Math.abs(b[0] - a[0]), 2) + Math.pow(Math.abs(b[1] - a[1]), 2))
    }
    addEdge(v1, v2) {
        // console.log('Connection',v1,'to',v2)
        this.adjList[v1].push(v2);
    }
    get conList() {
        const cl = [];
        let numEdges = 0;
        Object.keys(this.adjList).forEach(fromVert => {
            let pi = this.paths.find(q => q.id == fromVert);
            this.adjList[fromVert].forEach(toVert => {
                let pf = this.paths.find(q => q.id == toVert)
                cl.push({
                    from: fromVert.toString(),
                    to: toVert.toString(),
                    weight: this.dist(pi.center, pf.center),
                    visited: false,
                    used: false,
                    id: numEdges
                });
                numEdges++;
            })
        })
        return cl.sort((a, b) => a.weight - b.weight);
    }
    isConnected(cl, e) {
        const visited = [cl[0].from, cl[0].to];

        let limiter = 1000,
            currFroms = [...visited];

        while(limiter--){
            
        }
        // while(canContinue && limiter--){
        //     currFroms  = currFroms.map(cf=>{
        //         return cl.filter(q=>q.from==cf).map(c=>c.to);
        //     }).flat();
        //     console.log('Current Froms',currFroms,'used connections',cl,cl.length)
        //     if(currFroms.filter(cf=>visited.includes(cf)).length){
        //         //one of our outputs this round is already in our "visited" list, indicating a cycle!
        //         return false;
        //     }
        // }
        return true;
    }
    mst(cl) {
        cl[0].visited = true;
        cl[0].used = true;
        let currEdgeNum = 1,
            limiter = 10000;//just in case!
        while (limiter-- && currEdgeNum < cl.length && cl.filter(c => !c.visited)) {
            const currEdge = cl[currEdgeNum];
            currEdge.visited = true;//whether we add or not, set this to visited
            console.log('Examining edge #', currEdge.id, `(${currEdge.from} to ${currEdge.to})`)
            //if edge is not already connected (i.e., adding this would NOT create a cycle), mark it as "used"
            currEdge.used = !this.isConnected(cl.filter(q => !!q.used || q.id == currEdge.id), currEdge)
            currEdgeNum++;
        }
        return cl.filter(q => q.used);
        //at end, return all CL els that are used==true;
        // return [activeCons, visited];
    }
}